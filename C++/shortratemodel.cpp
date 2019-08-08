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
    TimeGrid times(maturity, timeSteps);
    boost::shared_ptr<OneFactorModel::ShortRateDynamics> CIRTree(CIRmodel_->dynamics());
    timeGrid_.resize(timeSteps + 1);
    std::copy(times.begin(), times.end(), timeGrid_.begin());
    Real r0 = CIRTree->shortRate(0.0, CIRTree->process()->x0());
    std::cout << r0 << std::endl;
    /*InterestRateTreeNode treeStructure(r0, numPaths);
    for(int i = 1; i < times.size(); i++){
        for(int j = 0; j < numPaths; j++){
            Real dW = 0.01;
            CIRTree.get()->process().get()->evolve();
        }
    } */
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