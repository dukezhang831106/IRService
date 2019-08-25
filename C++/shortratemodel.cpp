#include "shortratemodel.hpp"

HullWhiteModel::HullWhiteModel(std::string &modelDate, std::string &currency, std::string &exchange, std::string& interpolationType) : BaseModel(modelDate, currency, exchange, interpolationType) { 
    std::ifstream file("C++/data/configuration.json");
    json output;
    file >> output;
    for(auto & exchange : output[getCurrency()]){
        if (exchange["Exchange"] == getExchange()){
            // parse general information
            std::string holmask = exchange["Holidays"], settle = exchange["SettlementDays"], daycounter = exchange["DayCountConvention"];

            for (auto & instrument : exchange["VolatilityInstruments"].items()){
                std::string key(instrument.key());
                parse_fitting_bucket(key, instrument.value());
            }

            for (auto & instrument : exchange["VolatilityInstruments"].items()){
                std::string key(instrument.key());
                parse_vol_instruments(key, instrument.value());
            }
            buildModel();
        }
    }
    return;
};

void HullWhiteModel::parse_fitting_bucket(std::string& key, json& info){
    if (key == "CapFloorFittingGrid"){

    }
    else if (key == "SwaptionFittingGrid"){
        for (auto & vol_instrument : info.items()){
            std::string expiry = vol_instrument.value()["Expiry"], tenor = vol_instrument.value()["Tenor"];
            fittable_.insert("Swaption_" + expiry + "X" + tenor);
        }
    }
    return;
}

void HullWhiteModel::parse_vol_instruments(std::string& key, json& info){
    pqxx::connection conn("dbname=EikonInstrument user=postgres host=127.0.0.1 port=1111 password=123456");
    pqxx::work txn(conn);

    std::stringstream command;
    command << "SELECT * " << "FROM \"IRInstruments\" " << "WHERE \"QuoteDate\"=Date('" 
        << getModelDate() << "') AND \"Currency\"='USD' AND \"InstrumentType\"='" << key << "';";
    pqxx::result res = txn.exec(command.str());
    
    if (key == "Swaption"){
        std::string fix_dct = info["FixedLegDayCountConvention"], fix_bdc = info["FixedLegBusinessDayConvention"], fix_freq = info["FixedLegFrequency"];
        std::string flt_dct = info["FloatingLegDayCountConvention"], flt_bdc = info["FloatingLegBusinessDayConvention"], flt_freq = info["FloatingLegFrequency"];
        std::string settle = info["SettlementDays"], indices = info["Index"];
        
        bool eom = info["EOM"];
        int settlement_days = boost::lexical_cast<int>(settle.substr(0, settle.size() - 1));
        
        HWmodel_ = boost::shared_ptr<HullWhite>(new HullWhite(getDiscountTermStructure()));
        boost::shared_ptr<IborIndex> index = makeIndex(indices);
        boost::shared_ptr<PricingEngine> engine(new JamshidianSwaptionEngine(HWmodel_));

        for(int i = 0; i < res.size(); i++){
            pqxx::row record = res[i];
            GeneralInstrumentInformation instrument(record);
            std::string expiryTenor = "Swaption_" + instrument.expiry + "X" + instrument.tenor;
            if (fittable_.find(expiryTenor) != fittable_.end()){
                Period expiry = parse_period(instrument.expiry), tenor = parse_period(instrument.tenor);
                Handle<Quote> swaptionVol(boost::shared_ptr<Quote>(new SimpleQuote(instrument.quote/100.0)));

                boost::shared_ptr<CalibrationHelper> swaptionHelper(new SwaptionHelper(expiry, tenor, swaptionVol, index,
                                        index->tenor(), parse_daycounter(fix_dct),
                                        index->dayCounter(),
                                        getDiscountTermStructure(),
                                        CalibrationHelper::ImpliedVolError));

                swaptionHelper->setPricingEngine(engine);
                modelCalibrator_.push_back(swaptionHelper);
                CalibrationReport instr(key, instrument.expiry, instrument.tenor);
                calibrationReports_.push_back(instr);
            }
        }
        
    }
}

void HullWhiteModel::buildModel(){
    LevenbergMarquardt optimizationMethod(1.0e-8,1.0e-8,1.0e-8);
    EndCriteria endCriteria(10000, 100, 1e-6, 1e-8, 1e-8);

    //Optimize
    HWmodel_->calibrate(modelCalibrator_, optimizationMethod, endCriteria);
    EndCriteria::Type ecType = HWmodel_->endCriteria();

    Real calculated = 0.0;
    for (int i=0; i<modelCalibrator_.size(); ++i) {
        calibrationReports_[i].diff = modelCalibrator_[i]->calibrationError();
        calibrationReports_[i].calculated = modelCalibrator_[i]->modelValue();
        calibrationReports_[i].expected = modelCalibrator_[i]->marketValue();
        calculated += std::pow(calibrationReports_[i].diff, 2.0);
    }
    return;
}

void HullWhiteModel::generateMonteCarloPaths(int& numPaths, int& timeSteps, double& maturity){
    boost::shared_ptr<StochasticProcess1D> hwProcess(new HullWhiteProcess(getDiscountTermStructure(), HWmodel_.get()->a(), HWmodel_.get()->sigma()));
    // type definition for complex declaration
    typedef RandomSequenceGenerator<CLGaussianRng<MersenneTwisterUniformRng>> GSG;
    // create mersenne twister uniform random generator
    unsigned long seed = std::chrono::system_clock::now().time_since_epoch().count();;
    MersenneTwisterUniformRng generator(seed);
    // create gaussian generator by using central limit transformation method
    CLGaussianRng<MersenneTwisterUniformRng> gaussianGenerator(generator);
    GSG gaussianSequenceGenerator(timeSteps, gaussianGenerator);
    // create path generator using Hull-White process and gaussian sequence generator
    PathGenerator<GSG> pathGenerator(hwProcess, maturity, timeSteps, gaussianSequenceGenerator, false);
    TimeGrid times = pathGenerator.timeGrid();
    timeGrid_.resize(timeSteps + 1);
    std::copy(times.begin(), times.end(), timeGrid_.begin());
    paths_.resize(numPaths, timeSteps + 1);
    for(int i = 0; i < numPaths; i++){
        Sample<Path> IRpath = pathGenerator.next();
        for(int j = 0; j < IRpath.value.length(); j++){
            paths_(i, j) = IRpath.value.at(j);
        }
    }
    return;
}

void HullWhiteModel::buildProbabililtyTree(ProbabilityTree& ptree){
    int numLevels = ptree.R.size();
    boost::numeric::ublas::vector<double> V(numLevels), dr(numLevels+1); dr[0] = 0.0;
    for(int i = 0; i < numLevels; i++){
        V[i] = ptree.dT[i + 1]*std::pow(ptree.sigma, 2.0);
        dr[i+1] = ptree.sigma*std::sqrt(3*ptree.dT[i + 1]);
    }
    ptree.Branches[0] = boost::numeric::ublas::vector<int>(1, 0);
    
    for(int level = 0; level < numLevels; level++){
        // Calculate the rate values on this tree node
        ptree.R[level].resize(ptree.Branches[level].size());
        ptree.R[level] = dr[level]*ptree.Branches[level];
        if (level < numLevels - 1){
            boost::numeric::ublas::vector<double> M = -ptree.alpha*ptree.dT[level]*ptree.R[level];
            boost::numeric::ublas::vector<double> ConnectJ = (ptree.R[level] + M)/dr[level + 1];
            std::transform(ConnectJ.begin(), ConnectJ.end(), ConnectJ.begin(), [](double& x) {return std::round(x); });
            boost::numeric::ublas::vector<double> Al = (ptree.R[level] + M - dr[level + 1]*ConnectJ)/dr[level + 1];
            
            boost::numeric::ublas::matrix<double> prob(3, Al.size(), 0.0);
            boost::numeric::ublas::vector<double> pb0(boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<double>>(prob, 0));
            boost::numeric::ublas::vector<double> pb1(boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<double>>(prob, 1));
            boost::numeric::ublas::vector<double> pb2(boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<double>>(prob, 2));
 
            std::transform(Al.begin(), Al.end(), pb0.begin(), [&](double& x){ return (x*x + x)/2.0 + V[level]/std::pow(dr[level+1], 2.0)/2.0; });
            std::transform(Al.begin(), Al.end(), pb1.begin(), [&](double& x){ return 1.0 - std::pow(x, 2.0) - V[level]/std::pow(dr[level+1], 2.0); });
            std::transform(Al.begin(), Al.end(), pb2.begin(), [&](double& x){ return (x*x - x)/2.0 + V[level]/std::pow(dr[level+1], 2.0)/2.0; });
            boost::numeric::ublas::row(prob, 0) = pb0;
            boost::numeric::ublas::row(prob, 1) = pb1;
            boost::numeric::ublas::row(prob, 2) = pb2;

            ptree.Prob[level] = prob;

            std::set<double> conn;
            std::for_each(ConnectJ.begin(), ConnectJ.end(), [&](double x) { conn.insert(x); conn.insert(x+1); conn.insert(x-1); });
            std::vector<double> thisLevelJ(conn.size());
            std::partial_sort_copy(conn.begin(), conn.end(), thisLevelJ.begin(), thisLevelJ.end(), std::greater<int>());
            ptree.Branches[level+1].resize(thisLevelJ.size());
            std::copy(thisLevelJ.begin(), thisLevelJ.end(), ptree.Branches[level+1].begin());
            
            boost::numeric::ublas::matrix<int> connection(0, 0);
            for_each(ConnectJ.begin(), ConnectJ.end(), [&](double &x) { std::vector<double>::iterator it = std::find(thisLevelJ.begin(), thisLevelJ.end(), x); connection.resize(connection.size1()+1, 1); connection(connection.size1()-1, 0) = std::distance(thisLevelJ.begin(), it); });
            
            ptree.Connect[level] = connection;
        }   
    }
    ptree.dr = dr;
}


ProbabilityTree HullWhiteModel::getLattice(std::vector<Period>& maturities){
    Calendar calendar = getCalendar();
    std::string dc = getDayCounter();
    DayCounter daycounter = parse_daycounter(dc);
    Date valDate = getModelModelDate();
    std::vector<Date> volDates(maturities.size());
    std::transform(maturities.begin(), maturities.end(), volDates.begin(), [&](Period& p){ return calendar.advance(valDate, p); });
    Date alphaDate = volDates[volDates.size() - 1];

    double volCurve = HWmodel_->params()[1], alphaCurve = HWmodel_->params()[0];
    std::vector<double> timeFracs(volDates.size());
    std::transform(volDates.begin(), volDates.end(), timeFracs.begin(), [&](Date& d) { return daycounter.yearFraction(valDate, d); });
    std::vector<double> zRates = getCurveZeroRates(timeFracs), discountFactors = getCashDF(timeFracs);
    
    int numLevels = volDates.size();
    std::vector<double> tSpan(numLevels + 1, 0.0), dtSpan(numLevels + 1, 0.0); tSpan[0] = 0.0;
    std::copy(timeFracs.begin(), timeFracs.end(), tSpan.begin()+1);
    std::adjacent_difference(tSpan.begin(), tSpan.end(), dtSpan.begin());
    
    ProbabilityTree ptree(valDate, volDates, volCurve, alphaCurve, dtSpan);
    buildProbabililtyTree(ptree);
    
    boost::numeric::ublas::vector<double> shift(numLevels, 1.0); shift[0] = zRates[0];

    std::vector<boost::numeric::ublas::vector<double>> ADPTree = ptree.R;
    ADPTree[0] = boost::numeric::ublas::vector<double>(1, 1.0);
    ptree.R[0] = boost::numeric::ublas::vector<double>(1, shift[0]);
    for (int level = 1; level < numLevels; level++){
        boost::numeric::ublas::matrix<int> PrevConn(ptree.Connect[level - 1].size1(), 3, 1);
        boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<int>> prev(PrevConn, 1);
        boost::numeric::ublas::vector<int> pconn(boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<int>>(ptree.Connect[level - 1], 0));
        std::copy(pconn.begin(), pconn.end(), prev.begin());
        boost::numeric::ublas::column(PrevConn, 1) = prev;
        boost::numeric::ublas::column(PrevConn, 0) = boost::numeric::ublas::column(PrevConn, 1) - boost::numeric::ublas::column(PrevConn, 0);
        boost::numeric::ublas::column(PrevConn, 2) = boost::numeric::ublas::column(PrevConn, 1) + boost::numeric::ublas::column(PrevConn, 2);
        int jIdx = ADPTree[level].size();
        std::vector<double> temp(jIdx, 0);
        for(int iState = 0; iState < jIdx; iState++){
            std::vector<int> BackConn, direction;
            for(int row = 0; row < PrevConn.size1(); row++){
                for(int col = 0; col < PrevConn.size2(); col++){
                    if (PrevConn(row, col) == iState){
                        BackConn.push_back(row);
                        direction.push_back(col);
                    }
                }
            }
            for(int count = 0; count < BackConn.size(); count++){
                boost::numeric::ublas::matrix<double> prob = ptree.Prob[level - 1];
                boost::numeric::ublas::vector<double> R = ptree.R[level - 1], lastLevel = ADPTree[level - 1];
                double dt = ptree.dT[level];
                temp[iState] += prob(direction[count], BackConn[count]) * std::exp(-(ptree.R[level-1][BackConn[count]]) * dt) * lastLevel[BackConn[count]];
            }
            
        }
        std::copy(temp.begin(), temp.end(), ADPTree[level].begin());
        boost::numeric::ublas::vector<double> index(jIdx);
        std::iota(std::begin(index), std::end(index), 1.0);
        std::transform(index.begin(), index.end(), index.begin(), [&](double &x) { return std::floor(jIdx - x - jIdx/2); });

        boost::numeric::ublas::vector<double> expR(ptree.R[level].size());
        std::transform(ptree.R[level].begin(), ptree.R[level].end(), expR.begin(), [&](double& z) { return std::exp(-z*ptree.dT[level+1]); });
        shift[level] = std::log(boost::numeric::ublas::inner_prod(ADPTree[level], expR) - std::log(discountFactors[level]))/ptree.dT[level+1];
        //shift[level] = std::exp(-ptree.R[level]*ptree.dT[level+1]) - std::log(discountFactors[level])/ptree.dT[level+1];
        
        std::transform(ptree.R[level].begin(), ptree.R[level].end(), ptree.R[level].begin(), [&](double& x) { return x + shift[level]; });
    }

    // Convert to forward rates
    std::vector<boost::numeric::ublas::vector<double>> fwdRates = ptree.R;
    for (int level = 0; level < numLevels; level++){
        boost::numeric::ublas::vector<double> r = fwdRates[level];
        std::transform(r.begin(), r.end(), r.begin(), [](double& x) { return x; });
        std::transform(r.begin(), r.end(), fwdRates[level].begin(), [&](double& x) { return std::exp(x*ptree.dT[level+1]); });
        ptree.Fwd[level] = fwdRates[level];
    }
    
    return ptree;
}

CoxIngersollRossModel::CoxIngersollRossModel(std::string &modelDate, std::string &currency, std::string &exchange, std::string& interpolationType) : BaseModel(modelDate, currency, exchange, interpolationType) { 
    std::ifstream file("C++/data/configuration.json");
    json output;
    file >> output;
    for(auto & exchange : output[getCurrency()]){
        if (exchange["Exchange"] == getExchange()){
            // parse general information
            std::string holmask = exchange["Holidays"], settle = exchange["SettlementDays"], daycounter = exchange["DayCountConvention"];

            for (auto & instrument : exchange["VolatilityInstruments"].items()){
                std::string key(instrument.key());
                parse_fitting_bucket(key, instrument.value());
            }

            for (auto & instrument : exchange["VolatilityInstruments"].items()){
                std::string key(instrument.key());
                parse_vol_instruments(key, instrument.value());
            }
            buildModel();
        }
    }
    return;
};

void CoxIngersollRossModel::parse_fitting_bucket(std::string& key, json& info){
    if (key == "CapFloorFittingGrid"){

    }
    else if (key == "SwaptionFittingGrid"){
        for (auto & vol_instrument : info.items()){
            std::string expiry = vol_instrument.value()["Expiry"], tenor = vol_instrument.value()["Tenor"];
            fittable_.insert("Swaption_" + expiry + "X" + tenor);
        }
    }
    return;
}

void CoxIngersollRossModel::parse_vol_instruments(std::string& key, json& info){
    pqxx::connection conn("dbname=EikonInstrument user=postgres host=127.0.0.1 port=1111 password=123456");
    pqxx::work txn(conn);

    std::stringstream command;
    command << "SELECT * " << "FROM \"IRInstruments\" " << "WHERE \"QuoteDate\"=Date('" 
        << getModelDate() << "') AND \"Currency\"='USD' AND \"InstrumentType\"='" << key << "';";
    pqxx::result res = txn.exec(command.str());
    
    if (key == "Swaption"){
        std::string fix_dct = info["FixedLegDayCountConvention"], fix_bdc = info["FixedLegBusinessDayConvention"], fix_freq = info["FixedLegFrequency"];
        std::string flt_dct = info["FloatingLegDayCountConvention"], flt_bdc = info["FloatingLegBusinessDayConvention"], flt_freq = info["FloatingLegFrequency"];
        std::string settle = info["SettlementDays"], indices = info["Index"];
        
        bool eom = info["EOM"];
        int settlement_days = boost::lexical_cast<int>(settle.substr(0, settle.size() - 1));
        
        CIRmodel_ = boost::shared_ptr<ExtendedCoxIngersollRoss>(new ExtendedCoxIngersollRoss(getDiscountTermStructure()));
        boost::shared_ptr<IborIndex> index = makeIndex(indices);
        boost::shared_ptr<PricingEngine> engine(new JamshidianSwaptionEngine(CIRmodel_));

        for(int i = 0; i < res.size(); i++){
            pqxx::row record = res[i];
            GeneralInstrumentInformation instrument(record);
            std::string expiryTenor = "Swaption_" + instrument.expiry + "X" + instrument.tenor;
            if (fittable_.find(expiryTenor) != fittable_.end()){
                Period expiry = parse_period(instrument.expiry), tenor = parse_period(instrument.tenor);
                Handle<Quote> swaptionVol(boost::shared_ptr<Quote>(new SimpleQuote(instrument.quote/100.0)));

                boost::shared_ptr<CalibrationHelper> swaptionHelper(new SwaptionHelper(expiry, tenor, swaptionVol, index,
                                        index->tenor(), parse_daycounter(fix_dct),
                                        index->dayCounter(),
                                        getDiscountTermStructure(),
                                        CalibrationHelper::ImpliedVolError));

                swaptionHelper->setPricingEngine(engine);
                modelCalibrator_.push_back(swaptionHelper);
                CalibrationReport instr(key, instrument.expiry, instrument.tenor);
                calibrationReports_.push_back(instr);
            }
        }
        
    }
}

void CoxIngersollRossModel::buildModel(){
    LevenbergMarquardt optimizationMethod(1.0e-8,1.0e-8,1.0e-8);
    EndCriteria endCriteria(10000, 100, 1e-6, 1e-8, 1e-8);

    //Optimize
    CIRmodel_->calibrate(modelCalibrator_, optimizationMethod, endCriteria);
    EndCriteria::Type ecType = CIRmodel_->endCriteria();

    Real calculated = 0.0;
    for (int i=0; i<modelCalibrator_.size(); ++i) {
        calibrationReports_[i].diff = modelCalibrator_[i]->calibrationError();
        calibrationReports_[i].calculated = modelCalibrator_[i]->modelValue();
        calibrationReports_[i].expected = modelCalibrator_[i]->marketValue();
        calculated += std::pow(calibrationReports_[i].diff, 2.0);
    }
    return;
}

double CoxIngersollRossModel::findb(double& dT, double& X){
    double sqdt = std::sqrt(dT), xsqrtdt = X/std::sqrt(1.5*dT), bcom = X/sqdt;
    double base = xsqrtdt > 0 ? std::floor(xsqrtdt): -std::floor(-xsqrtdt), base2 = xsqrtdt+1 > 0 ? std::floor(xsqrtdt+1): -std::floor(-xsqrtdt+1);
    double be = bcom/base, bc = bcom/base2;
    return std::abs(bc - std::sqrt(1.5)) < std::abs(be - std::sqrt(1.5)) ? bc : be;
}

std::vector<double> CoxIngersollRossModel::solvedT(double& dT, double& X){
    int cnt = 0, Kmax = 100;
    double tol = 1e-10, prev0 = dT, curr0 = dT + tol, fprev0 = findb(prev0, X), fcurr0 = findb(curr0, X), prev1 = prev0, curr1 = curr0, fprev1 = fprev0, fcurr1 = fcurr0;
    std::vector<double> retval(2, 0.0);
    while(cnt < Kmax){
        double error0 = fcurr0 - std::sqrt(1.5), error1 = fcurr1 - 1.0;
        if (std::abs(error0) < tol && std::abs(error1) < tol){
            retval[0] = curr0;
            retval[1] = curr1;
            return retval;
        }
        else{
            if (std::abs(error0) > tol){
                double temp = curr0 - (curr0 - prev0) * error0/(fcurr0 - fprev0);
                prev0 = curr0;
                fprev0 = fcurr0;
                curr0 = temp;
                fcurr0 = findb(curr0, X);
            }
            if (std::abs(error1) > tol){
                double temp = curr1 - (curr1 - prev1) * error1/(fcurr1 - fprev1);
                prev1 = curr1;
                fprev1 = fcurr1;
                curr1 = temp;
                fcurr1 = findb(curr1, X);
            }
        }
    }
    retval[0] = std::numeric_limits<double>::infinity();
    retval[1] = std::numeric_limits<double>::infinity();
    return retval;
}

// remove elements from a vector if value is approximately near (up to 9 decimals)
void CoxIngersollRossModel::remove_duplicates(std::vector<double>& vec){
    int rounder = 1e8;
    std::map<double, double> mapping;
    std::for_each(vec.begin(), vec.end(), [&](double& x) {
        double key = std::round(rounder*x)/rounder;
        if (mapping.find(key) == mapping.end()){
            mapping[key] = x;
        }
            
    });
    vec.resize(mapping.size());
    int cnt = 0;
    for(std::map<double, double>::iterator it = mapping.begin(); it != mapping.end(); it++){
        vec[cnt] = it->first; cnt++;
    }
    std::sort(vec.begin(), vec.end(), std::greater<double>());
}

ProbabilityTree CoxIngersollRossModel::getLattice(std::vector<Period>& maturities){
    Calendar calendar = getCalendar();
    std::string dct = getDayCounter();
    DayCounter daycounter = Actual365NoLeap();
    Date valDate = getModelModelDate();
    std::vector<Date> volDates(maturities.size() + 1);
    volDates[0] = valDate;
    std::transform(maturities.begin(), maturities.end(), volDates.begin() + 1, [&](Period& p){ return calendar.advance(valDate, p); });
    Date alphaDate = volDates[volDates.size() - 1];

    double thetaCurve = CIRmodel_->params()[0], alphaCurve = CIRmodel_->params()[1], volCurve = CIRmodel_->params()[2];
    std::vector<double> timeFracs(volDates.size());
    std::transform(volDates.begin(), volDates.end(), timeFracs.begin(), [&](Date& d) { return daycounter.yearFraction(valDate, d); });
    
    std::vector<double> zRates = getCurveZeroRates(timeFracs);
    std::vector<double> discountFactors(zRates.size()), nominalBond(zRates.size());
    int timeSteps = -1;
    std::transform(zRates.begin(), zRates.end(), discountFactors.begin(), [&](double& r) { timeSteps++; return std::pow(1+r, -timeFracs[timeSteps]); });

    double dT = timeFracs[timeFracs.size() - 1] / timeFracs.size();
    
    // linear interpolation the first grid
    LinearInterpolation Interp (timeFracs.begin(),timeFracs.end(),zRates.begin());

    double refRate = Interp(dT);
    timeSteps = 0;
    // Find Prices of nominal bond
    std::for_each(timeFracs.begin(), timeFracs.end(), [&](double& t) { timeSteps++; double rate = Interp(timeSteps*dT); nominalBond[timeSteps-1] = std::pow(1 + rate, -timeSteps*dT); });
    
    int numLevels = volDates.size();
    std::vector<double> tSpan(numLevels + 1, 0.0), dtSpan(numLevels + 1, 0.0); tSpan[0] = 0.0;
    std::copy(timeFracs.begin(), timeFracs.end(), tSpan.begin()+1);
    std::adjacent_difference(tSpan.begin(), tSpan.end(), dtSpan.begin());
    
    ProbabilityTree ptree(valDate, volDates, volCurve, thetaCurve, alphaCurve, dtSpan);

    double Xprev = 2*sqrt(refRate)/ptree.sigma, sqdt = std::sqrt(dT), b = findb(dT, Xprev);
    std::vector<double> vecXprev(1, Xprev);

    if (b > std::sqrt(2.0) || b < 1.0){
        // find dT for which b will make sense
        std::vector<double> dt = solvedT(dT, Xprev);

    }

    double sigsq = ptree.sigma * ptree.sigma, sigsq4 = sigsq / 4.0;

    ptree.R[0] = boost::numeric::ublas::vector<double>(1, refRate);
    // jump if necessary
    double Xju = 2*std::sqrt(ptree.alpha*ptree.theta*dT)/ptree.sigma;


    for (int level = 1; level < numLevels; level++){
        std::vector<double> Xu(vecXprev.size()), Xm(vecXprev.size()), Xd(vecXprev.size());
        std::vector<int> nomask;       
        for (int index = 0; index < vecXprev.size(); index++){
            double muxt = (0.5 * ptree.alpha * (4.0 * ptree.theta / sigsq - vecXprev[index] * vecXprev[index]) - 0.5) / vecXprev[index];
            double J = std::floor(muxt * sqdt / b + 1.0 / b / b);
            if (vecXprev[index] > DBL_EPSILON){    
                Xu[index] = vecXprev[index] + b * (J + 1) * sqdt;
                Xm[index] = vecXprev[index] + b * J * sqdt;
                Xd[index] = std::max(vecXprev[index] + b * (J - 1) * sqdt, 0.0);
                Xd[index] = Xd[index] < 100000 * DBL_EPSILON ? 0.0 : Xd[index];
            }
            else{
                Xu[index] = 0.0;
                Xm[index] = 0.0;
                Xd[index] = 0.0;
                nomask.push_back(index);
            }
            ptree.Prob[level - 1].resize(3, vecXprev.size());
            ptree.Prob[level - 1](0, index) = 0.5/b/b - J/2 + sqdt*muxt/2/b;
            ptree.Prob[level - 1](1, index) = 1 - 1/b/b;
            ptree.Prob[level - 1](2, index) = 0.5/b/b + J/2 - sqdt*muxt/2/b;
        }
        remove_duplicates(Xu);
        remove_duplicates(Xm);
        remove_duplicates(Xd);

        std::vector<double> vecX(Xu);
        vecX.insert(vecX.end(), Xm.begin(), Xm.end());
        vecX.insert(vecX.end(), Xd.begin(), Xd.end());

        remove_duplicates(vecX);

        // update rates

        ptree.R[level].resize(vecX.size());
        std::transform(vecX.begin(), vecX.end(), ptree.R[level].begin(), [&](double& x) { return sigsq4*x*x; });
        int idxJump = -1;
        if (std::abs(ptree.R[level-1][ptree.R[level-1].size()-1]) < DBL_EPSILON){
            std::vector<double>::reverse_iterator low = std::lower_bound(Xu.rbegin(), Xu.rend(), Xju, std::less<double>());
            idxJump = std::distance(low, Xu.rend()) - 1;
        }


        // update Connect
        boost::numeric::ublas::matrix<int> conn(3, Xu.size());
        
        for(int i = 0; i < Xu.size(); i++){
            std::vector<double>::iterator it = std::find(vecX.begin(), vecX.end(), Xu[i]); 
            conn(0, i) = std::distance(vecX.begin(), it);
            it = std::find(vecX.begin(), vecX.end(), Xm[i]); 
            conn(1, i) = std::distance(vecX.begin(), it);
            it = std::find(vecX.begin(), vecX.end(), Xd[i]); 
            conn(2, i) = std::distance(vecX.begin(), it);
            if (idxJump >= 0){
                std::for_each(nomask.begin(), nomask.end(), [&](int& pos) { conn(0, pos) = idxJump;});
            }
            
        }
        ptree.Connect[level-1] = conn;

        if (idxJump >= 0){
            boost::numeric::ublas::vector<double> pu(conn.size2()), temp_u(conn.size2()), temp_d(conn.size2());
            boost::numeric::ublas::vector<int> puconn(boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<int>>(conn, 0));
            boost::numeric::ublas::vector<int> pdconn(boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<int>>(conn, 2));
            std::transform(ptree.R[level-1].begin(), ptree.R[level-1].end(), pu.begin(), [&](double& r) { return r + dT*ptree.alpha*(ptree.theta - r); });
            std::transform(puconn.begin(), puconn.end(), temp_u.begin(), [&](int& idx) { return ptree.R[level][idx]; });
            std::transform(pdconn.begin(), pdconn.end(), temp_d.begin(), [&](int& idx) { return ptree.R[level][idx]; });
            pu = pu - temp_d;
            temp_u = temp_u - temp_d;
            std::transform(pu.begin(), pu.end(), temp_u.begin(), pu.begin(), std::divides<double>());
            std::for_each(nomask.begin(), nomask.end(), [&](int& pos) { ptree.Prob[level-1](0, pos) = pu(pos); ptree.Prob[level-1](1, pos) = 0.0; ptree.Prob[level-1](2, pos) = 1.0 - pu(pos);});
        }
        
        // update probabilities
        vecXprev = vecX;
    }
    
    // Now calculate forward rates
    ptree.Fwd = ptree.R;
    for(int level = 0; level < ptree.R.size(); level++){
        std::transform(ptree.R[level].begin(), ptree.R[level].end(), ptree.Fwd[level].begin(), [&](double& r) { return std::exp(r*dT); });
    }

    // The trees we have are for the initial rate (first time Step). We have to now
    // adjust the tree to the rest of the term structure. We do that by
    // calculating the price of pseudo-bonds with face=1 maturing at each tree level

    std::vector<boost::numeric::ublas::vector<double>> Pt(ptree.Fwd.size() + 1);
    std::vector<boost::numeric::ublas::matrix<double>> ProbT(Pt.size());
    std::copy(ptree.Fwd.begin(), ptree.Fwd.end(), Pt.begin());
    std::copy(ptree.Fwd.end() - 1, ptree.Fwd.end(), Pt.end() - 1);
    
    std::transform(Pt.begin(), Pt.end(), ProbT.begin(), [&](boost::numeric::ublas::vector<double>& row) { return boost::numeric::ublas::matrix<double>(numLevels, row.size(), 0.0); });
    
    boost::numeric::ublas::matrix<double> Face(numLevels , numLevels + 1, 0.0);
    boost::numeric::ublas::subrange(Face, 0, numLevels, 1, numLevels + 1) = boost::numeric::ublas::identity_matrix<double>(numLevels);
   
    int NPLevels = Pt.size(), NBranches = Pt[1].size();
    for(int level = NPLevels - 1; level >= 0; level--){
        if (level < NPLevels - 1){
            for(int i = 0; i < ProbT[level].size1(); i++){
                for(int j = 0; j < ProbT[level].size2(); j++)
                    ProbT[level](i, j) = ProbT[level](i, j)/ptree.Fwd[level][j];
            }
        }
        
        ProbT[level] += boost::numeric::ublas::prod(boost::numeric::ublas::subrange(Face, 0, numLevels, level, level + 1), boost::numeric::ublas::matrix<double>(1, Pt[level].size(), 1.0)); 

        if (level == NPLevels - 1){
            ProbT[level - 1] = ProbT[level];
        }
        else if (level >= 1){
            for(int branch = 0; branch < NBranches; branch++){
                boost::numeric::ublas::vector<int> connectTo(boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<int>>(ptree.Connect[level - 1], branch));
                for(int i = 0; i < ProbT[level - 1].size1(); i++){
                   for(int j = 0; j < ProbT[level - 1].size2(); j++){
                        ProbT[level - 1](i, j) += ptree.Prob[level - 1](branch, j) * ProbT[level](i, connectTo[j]);
                   }
                }
            }
        }
    }

    boost::numeric::ublas::vector<double> Pdt(discountFactors.size()), delta(Pdt.size());
    std::transform(ProbT[0].begin1(), ProbT[0].end1(), nominalBond.begin(), Pdt.begin(), std::divides<double>());
    std::transform(Pdt.begin(), Pdt.end(), Pdt.begin(), [&](double& z) { return std::log(z)/dT; });
    
    std::adjacent_difference(Pdt.begin(), Pdt.end(), delta.begin());

    for(int level = 0; level < numLevels; level++){
        std::transform(ptree.R[level].begin(), ptree.R[level].end(), ptree.R[level].begin(), [&](double& r) { return r + delta[level]; });
        std::transform(ptree.R[level].begin(), ptree.R[level].end(), ptree.Fwd[level].begin(), [&](double& r) { return std::exp(r*dT); });
    }

    
    return ptree;
}


BlackKarasinskiModel::BlackKarasinskiModel(std::string &modelDate, std::string &currency, std::string &exchange, std::string& interpolationType) : BaseModel(modelDate, currency, exchange, interpolationType) { 
    std::ifstream file("C++/data/configuration.json");
    json output;
    file >> output;
    for(auto & exchange : output[getCurrency()]){
        if (exchange["Exchange"] == getExchange()){
            // parse general information
            std::string holmask = exchange["Holidays"], settle = exchange["SettlementDays"], daycounter = exchange["DayCountConvention"];

            for (auto & instrument : exchange["VolatilityInstruments"].items()){
                std::string key(instrument.key());
                parse_fitting_bucket(key, instrument.value());
            }

            for (auto & instrument : exchange["VolatilityInstruments"].items()){
                std::string key(instrument.key());
                parse_vol_instruments(key, instrument.value());
            }
            buildModel();
        }
    }
    return;
};

void BlackKarasinskiModel::parse_fitting_bucket(std::string& key, json& info){
    if (key == "CapFloorFittingGrid"){

    }
    else if (key == "SwaptionFittingGrid"){
        for (auto & vol_instrument : info.items()){
            std::string expiry = vol_instrument.value()["Expiry"], tenor = vol_instrument.value()["Tenor"];
            fittable_.insert("Swaption_" + expiry + "X" + tenor);
        }
    }
    return;
}

void BlackKarasinskiModel::parse_vol_instruments(std::string& key, json& info){
    pqxx::connection conn("dbname=EikonInstrument user=postgres host=127.0.0.1 port=1111 password=123456");
    pqxx::work txn(conn);

    std::stringstream command;
    command << "SELECT * " << "FROM \"IRInstruments\" " << "WHERE \"QuoteDate\"=Date('" 
        << getModelDate() << "') AND \"Currency\"='USD' AND \"InstrumentType\"='" << key << "';";
    pqxx::result res = txn.exec(command.str());
    
    if (key == "Swaption"){
        std::string fix_dct = info["FixedLegDayCountConvention"], fix_bdc = info["FixedLegBusinessDayConvention"], fix_freq = info["FixedLegFrequency"];
        std::string flt_dct = info["FloatingLegDayCountConvention"], flt_bdc = info["FloatingLegBusinessDayConvention"], flt_freq = info["FloatingLegFrequency"];
        std::string settle = info["SettlementDays"], indices = info["Index"];
        
        bool eom = info["EOM"];
        int settlement_days = boost::lexical_cast<int>(settle.substr(0, settle.size() - 1));
        
        BKmodel_ = boost::shared_ptr<BlackKarasinski>(new BlackKarasinski(getDiscountTermStructure()));
        boost::shared_ptr<IborIndex> index = makeIndex(indices);
        Size grid = 100;
        boost::shared_ptr<PricingEngine> engine(new TreeSwaptionEngine(BKmodel_, grid, getDiscountTermStructure()));

        for(int i = 0; i < res.size(); i++){
            pqxx::row record = res[i];
            GeneralInstrumentInformation instrument(record);
            std::string expiryTenor = "Swaption_" + instrument.expiry + "X" + instrument.tenor;
            if (fittable_.find(expiryTenor) != fittable_.end()){
                Period expiry = parse_period(instrument.expiry), tenor = parse_period(instrument.tenor);
                Handle<Quote> swaptionVol(boost::shared_ptr<Quote>(new SimpleQuote(instrument.quote/100.0)));

                boost::shared_ptr<CalibrationHelper> swaptionHelper(new SwaptionHelper(expiry, tenor, swaptionVol, index,
                                        index->tenor(), parse_daycounter(fix_dct),
                                        index->dayCounter(),
                                        getDiscountTermStructure(),
                                        CalibrationHelper::ImpliedVolError));

                swaptionHelper->setPricingEngine(engine);
                modelCalibrator_.push_back(swaptionHelper);
                CalibrationReport instr(key, instrument.expiry, instrument.tenor);
                calibrationReports_.push_back(instr);
            }
        }
        
    }
}

void BlackKarasinskiModel::buildModel(){
    LevenbergMarquardt optimizationMethod(1.0e-8,1.0e-8,1.0e-8);
    EndCriteria endCriteria(10000, 100, 1e-6, 1e-8, 1e-8);

    //Optimize
    BKmodel_->calibrate(modelCalibrator_, optimizationMethod, endCriteria);
    EndCriteria::Type ecType = BKmodel_->endCriteria();

    Real calculated = 0.0;
    for (int i=0; i<modelCalibrator_.size(); ++i) {
        calibrationReports_[i].diff = modelCalibrator_[i]->calibrationError();
        calibrationReports_[i].calculated = modelCalibrator_[i]->modelValue();
        calibrationReports_[i].expected = modelCalibrator_[i]->marketValue();
        calculated += std::pow(calibrationReports_[i].diff, 2.0);
    }
    return;
}

void BlackKarasinskiModel::buildProbabililtyTree(ProbabilityTree& ptree){
    int numLevels = ptree.R.size();
    boost::numeric::ublas::vector<double> V(numLevels), dr(numLevels+1); dr[0] = 0.0;
    for(int i = 0; i < numLevels; i++){
        V[i] = ptree.dT[i + 1]*std::pow(ptree.sigma, 2.0);
        dr[i+1] = ptree.sigma*std::sqrt(3*ptree.dT[i + 1]);
    }
    ptree.Branches[0] = boost::numeric::ublas::vector<int>(1, 0);
    
    for(int level = 0; level < numLevels; level++){
        // Calculate the rate values on this tree node
        ptree.R[level].resize(ptree.Branches[level].size());
        ptree.R[level] = dr[level]*ptree.Branches[level];
        if (level < numLevels - 1){
            boost::numeric::ublas::vector<double> M = -ptree.alpha*ptree.dT[level]*ptree.R[level];
            boost::numeric::ublas::vector<double> ConnectJ = (ptree.R[level] + M)/dr[level + 1];
            std::transform(ConnectJ.begin(), ConnectJ.end(), ConnectJ.begin(), [](double& x) {return std::round(x); });
            boost::numeric::ublas::vector<double> Al = (ptree.R[level] + M - dr[level + 1]*ConnectJ)/dr[level + 1];
            
            boost::numeric::ublas::matrix<double> prob(3, Al.size(), 0.0);
            boost::numeric::ublas::vector<double> pb0(boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<double>>(prob, 0));
            boost::numeric::ublas::vector<double> pb1(boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<double>>(prob, 1));
            boost::numeric::ublas::vector<double> pb2(boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<double>>(prob, 2));
 
            std::transform(Al.begin(), Al.end(), pb0.begin(), [&](double& x){ return (x*x + x)/2.0 + V[level]/std::pow(dr[level+1], 2.0)/2.0; });
            std::transform(Al.begin(), Al.end(), pb1.begin(), [&](double& x){ return 1.0 - std::pow(x, 2.0) - V[level]/std::pow(dr[level+1], 2.0); });
            std::transform(Al.begin(), Al.end(), pb2.begin(), [&](double& x){ return (x*x - x)/2.0 + V[level]/std::pow(dr[level+1], 2.0)/2.0; });
            boost::numeric::ublas::row(prob, 0) = pb0;
            boost::numeric::ublas::row(prob, 1) = pb1;
            boost::numeric::ublas::row(prob, 2) = pb2;

            ptree.Prob[level] = prob;

            std::set<double> conn;
            std::for_each(ConnectJ.begin(), ConnectJ.end(), [&](double x) { conn.insert(x); conn.insert(x+1); conn.insert(x-1); });
            std::vector<double> thisLevelJ(conn.size());
            std::partial_sort_copy(conn.begin(), conn.end(), thisLevelJ.begin(), thisLevelJ.end(), std::greater<int>());
            ptree.Branches[level+1].resize(thisLevelJ.size());
            std::copy(thisLevelJ.begin(), thisLevelJ.end(), ptree.Branches[level+1].begin());
            
            boost::numeric::ublas::matrix<int> connection(0, 0);
            for_each(ConnectJ.begin(), ConnectJ.end(), [&](double &x) { std::vector<double>::iterator it = std::find(thisLevelJ.begin(), thisLevelJ.end(), x); connection.resize(connection.size1()+1, 1); connection(connection.size1()-1, 0) = std::distance(thisLevelJ.begin(), it); });
            
            ptree.Connect[level] = connection;
        }   
    }
    ptree.dr = dr;
}

// solve for BK tree shift by Newton's method
double BlackKarasinskiModel::fsolve(double& init, double& targetPrice , boost::numeric::ublas::vector<double>& branches, boost::numeric::ublas::vector<double>& treeLevel, double& rateDelta, double& timeDelta){
    double EPS = 1e-10, prev = init, error = 1, curr;
    int Kmax = 100, cnt = 0;
    while(error > EPS){
        boost::numeric::ublas::vector<double> temp(branches.size());
        std::transform(branches.begin(), branches.end(), temp.begin(), [&](double& j) {return std::exp(-timeDelta*std::exp(prev + j*rateDelta)); });
        double fprev = boost::numeric::ublas::inner_prod(temp, treeLevel) - targetPrice;
        std::transform(branches.begin(), branches.end(), temp.begin(), [&](double& j) {return -timeDelta*std::exp(j*rateDelta + prev -timeDelta*std::exp(prev + j*rateDelta)); });
        double dfprev = boost::numeric::ublas::inner_prod(temp, treeLevel);
        curr = prev - fprev/dfprev;
        error = std::abs(fprev/dfprev);
        if (cnt > Kmax){
            return std::numeric_limits<double>::infinity();
        }
        cnt++;
        prev = curr;
    }
    return curr;
}

ProbabilityTree BlackKarasinskiModel::getLattice(std::vector<Period>& maturities){
    Calendar calendar = getCalendar();
    std::string dc = getDayCounter();
    DayCounter daycounter = parse_daycounter(dc);
    Date valDate = getModelModelDate();
    std::vector<Date> volDates(maturities.size());
    std::transform(maturities.begin(), maturities.end(), volDates.begin(), [&](Period& p){ return calendar.advance(valDate, p); });
    Date alphaDate = volDates[volDates.size() - 1];

    double volCurve = BKmodel_->params()[1], alphaCurve = BKmodel_->params()[0];
    std::vector<double> timeFracs(volDates.size());
    std::transform(volDates.begin(), volDates.end(), timeFracs.begin(), [&](Date& d) { return daycounter.yearFraction(valDate, d); });
    std::vector<double> zRates = getCurveZeroRates(timeFracs), discountFactors = getCashDF(timeFracs);
    
    int numLevels = volDates.size();
    std::vector<double> tSpan(numLevels + 1, 0.0), dtSpan(numLevels + 1, 0.0); tSpan[0] = 0.0;
    std::copy(timeFracs.begin(), timeFracs.end(), tSpan.begin()+1);
    std::adjacent_difference(tSpan.begin(), tSpan.end(), dtSpan.begin());
    
    ProbabilityTree ptree(valDate, volDates, volCurve, alphaCurve, dtSpan);
    buildProbabililtyTree(ptree);
    
    boost::numeric::ublas::vector<double> shift(numLevels, 1.0); shift[0] = std::log(zRates[0]);

    std::vector<boost::numeric::ublas::vector<double>> ADPTree = ptree.R;
    ADPTree[0] = boost::numeric::ublas::vector<double>(1, 1.0);
    ptree.R[0] = boost::numeric::ublas::vector<double>(1, shift[0]);
    for (int level = 1; level < numLevels; level++){
        boost::numeric::ublas::matrix<int> PrevConn(ptree.Connect[level - 1].size1(), 3, 1);
        boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<int>> prev(PrevConn, 1);
        boost::numeric::ublas::vector<int> pconn(boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<int>>(ptree.Connect[level - 1], 0));
        std::copy(pconn.begin(), pconn.end(), prev.begin());
        boost::numeric::ublas::column(PrevConn, 1) = prev;
        boost::numeric::ublas::column(PrevConn, 0) = boost::numeric::ublas::column(PrevConn, 1) - boost::numeric::ublas::column(PrevConn, 0);
        boost::numeric::ublas::column(PrevConn, 2) = boost::numeric::ublas::column(PrevConn, 1) + boost::numeric::ublas::column(PrevConn, 2);
        int jIdx = ADPTree[level].size();
        std::vector<double> temp(jIdx, 0);
        for(int iState = 0; iState < jIdx; iState++){
            std::vector<int> BackConn, direction;
            for(int row = 0; row < PrevConn.size1(); row++){
                for(int col = 0; col < PrevConn.size2(); col++){
                    if (PrevConn(row, col) == iState){
                        BackConn.push_back(row);
                        direction.push_back(col);
                    }
                }
            }
            for(int count = 0; count < BackConn.size(); count++){
                boost::numeric::ublas::matrix<double> prob = ptree.Prob[level - 1];
                boost::numeric::ublas::vector<double> R = ptree.R[level - 1], lastLevel = ADPTree[level - 1];
                double dt = ptree.dT[level];
                temp[iState] += prob(direction[count], BackConn[count]) * std::exp(-std::exp(ptree.R[level-1][BackConn[count]]) * dt) * lastLevel[BackConn[count]];
            }
            
        }
        std::copy(temp.begin(), temp.end(), ADPTree[level].begin());
        boost::numeric::ublas::vector<double> index(jIdx);
        std::iota(std::begin(index), std::end(index), 1.0);
        std::transform(index.begin(), index.end(), index.begin(), [&](double &x) { return std::floor(jIdx - x - jIdx/2); });

        shift[level] = fsolve(shift[0], discountFactors[level], index, ADPTree[level], ptree.dr[level], ptree.dT[level+1]);
        
        std::transform(ptree.R[level].begin(), ptree.R[level].end(), ptree.R[level].begin(), [&](double& x) { return x + shift[level]; });
    }

    // Convert to forward rates
    std::vector<boost::numeric::ublas::vector<double>> fwdRates = ptree.R;
    for (int level = 0; level < numLevels; level++){
        boost::numeric::ublas::vector<double> r = fwdRates[level];
        std::transform(r.begin(), r.end(), r.begin(), [](double& x) { return std::exp(x); });
        std::transform(r.begin(), r.end(), fwdRates[level].begin(), [&](double& x) { return std::exp(x*ptree.dT[level+1]); });
        ptree.Fwd[level] = fwdRates[level];
    }
    
    return ptree;
}