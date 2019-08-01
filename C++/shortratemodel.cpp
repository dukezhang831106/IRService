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
            pqxx::result::tuple record = res[i];
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
                std::cout << expiryTenor << " inserted." << std::endl;
                modelCalibrator_.push_back(swaptionHelper);
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

    std::cout << HWmodel_.get()->params() << std::endl;
    return;
}

ExtendedCIRModel::ExtendedCIRModel(std::string &modelDate, std::string &currency, std::string &exchange, std::string& interpolationType) : BaseModel(modelDate, currency, exchange, interpolationType) { 
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

void ExtendedCIRModel::parse_fitting_bucket(std::string& key, json& info){
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

void ExtendedCIRModel::parse_vol_instruments(std::string& key, json& info){
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
        
        eCIRmodel_ = boost::shared_ptr<ExtendedCoxIngersollRoss>(new ExtendedCoxIngersollRoss(getDiscountTermStructure()));
        boost::shared_ptr<IborIndex> index = makeIndex(indices);
        boost::shared_ptr<PricingEngine> engine(new JamshidianSwaptionEngine(eCIRmodel_));

        for(int i = 0; i < res.size(); i++){
            pqxx::result::tuple record = res[i];
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
                std::cout << expiryTenor << " inserted." << std::endl;
                modelCalibrator_.push_back(swaptionHelper);
            }
        }
        
    }
}

void ExtendedCIRModel::buildModel(){
    LevenbergMarquardt optimizationMethod(1.0e-8,1.0e-8,1.0e-8);
    EndCriteria endCriteria(10000, 100, 1e-6, 1e-8, 1e-8);

    //Optimize
    eCIRmodel_->calibrate(modelCalibrator_, optimizationMethod, endCriteria);
    EndCriteria::Type ecType = eCIRmodel_->endCriteria();

    std::cout << eCIRmodel_.get()->params() << std::endl;
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
        /*std::vector<Time> modelTimes;
        Date d0 = getSettlementDate();
        for(auto& v : getCurveCalibrator()){
            Date d = v.get()->maturityDate();
            modelTimes.push_back(getDiscountTermStructure().currentLink().get()->dayCounter().yearFraction(d0, d));
        }
        TimeGrid grid(modelTimes.begin(), modelTimes.end());*/
        Size grid = 100;
        boost::shared_ptr<PricingEngine> engine(new TreeSwaptionEngine(BKmodel_, grid, getDiscountTermStructure()));

        for(int i = 0; i < res.size(); i++){
            pqxx::result::tuple record = res[i];
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
                std::cout << expiryTenor << " inserted." << std::endl;
                modelCalibrator_.push_back(swaptionHelper);
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

    std::cout << BKmodel_.get()->params() << std::endl;
    return;
}