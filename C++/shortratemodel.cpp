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
    pqxx::connection conn("dbname=EikonInstrument user=postgres host=127.0.0.1 port=5432 password=123456");
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
    //std::cout << timeGrid_ << std::endl;
    paths_.resize(numPaths, timeSteps + 1);
    for(int i = 0; i < numPaths; i++){
        Sample<Path> IRpath = pathGenerator.next();
        for(int j = 0; j < IRpath.value.length(); j++){
            paths_(i, j) = IRpath.value.at(j);
        }
    }
    //std::cout << paths_ << std::endl;
    return;
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
    pqxx::connection conn("dbname=EikonInstrument user=postgres host=127.0.0.1 port=5432 password=123456");
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

void CoxIngersollRossModel::generateMonteCarloPaths(int& numPaths, int& timeSteps, double& maturity){
    TimeGrid grid(maturity, timeSteps);
    timeGrid_.resize(timeSteps + 1);
    std::copy(grid.begin(), grid.end(), timeGrid_.begin());




    TermStructureFittingParameter phi(getDiscountTermStructure());
    boost::shared_ptr<Lattice> tree = CIRmodel_->tree(grid);
    for(int i = 0; i < grid.size(); i++){
        Time t = grid.at(i);
        std::cout << tree.get()->grid(t) << std::endl;
    }   
    return;
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
    pqxx::connection conn("dbname=EikonInstrument user=postgres host=127.0.0.1 port=5432 password=123456");
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
            
            std::vector<int> connection;
            for_each(ConnectJ.begin(), ConnectJ.end(), [&](double &x) { std::vector<double>::iterator it = std::find(thisLevelJ.begin(), thisLevelJ.end(), x); connection.push_back(std::distance(thisLevelJ.begin(), it)); });
            
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
        boost::numeric::ublas::matrix<int> PrevConn(ptree.Connect[level - 1].size(), 3, 1);
        boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<int>> prev(PrevConn, 1);
        std::copy(ptree.Connect[level - 1].begin(), ptree.Connect[level - 1].end(), prev.begin());
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
        ptree.R[level] = fwdRates[level];
    }
    
    return ptree;
}