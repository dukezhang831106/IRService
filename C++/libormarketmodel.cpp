#include "libormarketmodel.hpp"

LiborMarketModel::LiborMarketModel(std::string &modelDate, std::string &currency, std::string &exchange, std::string& interpolationType) : BaseModel(modelDate, currency, exchange, interpolationType) { 
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


void LiborMarketModel::parse_fitting_bucket(std::string& key, json& info){
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

void LiborMarketModel::parse_vol_instruments(std::string& key, json& info){
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
        boost::shared_ptr<IborIndex> index = makeIndex(indices);
        bool eom = info["EOM"];
        int settlement_days = boost::lexical_cast<int>(settle.substr(0, settle.size() - 1));
            
        const Size size = 120;
        boost::shared_ptr<LiborForwardModelProcess> process(new LiborForwardModelProcess(size, index));
        // set-up the model
        boost::shared_ptr<LmVolatilityModel> volModel(new LmExtLinearExponentialVolModel(process->fixingTimes(), 0.5,0.6,0.1,0.1));
        boost::shared_ptr<LmCorrelationModel> corrModel(new LmLinearExponentialCorrelationModel(size, 0.5, 0.8));
        LMMmodel_ = boost::shared_ptr<LiborForwardModel>(new LiborForwardModel(process, volModel, corrModel));

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

                swaptionHelper->setPricingEngine(
                        boost::shared_ptr<PricingEngine>(
                                    new LfmSwaptionEngine(LMMmodel_, getDiscountTermStructure())));
                    std::cout << expiryTenor << " inserted." << std::endl;
                    modelCalibrator_.push_back(swaptionHelper); 
            }
        }
        
    }
}

void LiborMarketModel::buildModel(){
    LevenbergMarquardt om(1e-6, 1e-6, 1e-6);
    LMMmodel_->calibrate(modelCalibrator_, om, EndCriteria(2000, 100, 1e-6, 1e-6, 1e-6));
    std::cout << "Model Constructed..." << std::endl;
    // measure the calibration error
        Real calculated = 0.0;
        for (int i=0; i<modelCalibrator_.size(); ++i) {
            Real diff = modelCalibrator_[i]->calibrationError();
            //Real mkt = modelCalibrator_[i]->marketValue();
            //Real mdl = modelCalibrator_[i]->modelValue();
            std::cout << "Diff: " << diff << std::endl;
            //std::cout << "Market: " << mkt << ", Model: " << mdl << ", Diff: " << diff << std::endl;
            calculated += diff*diff;
        } 
}