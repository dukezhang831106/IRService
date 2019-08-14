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
        index_ = makeIndex(indices);
        bool eom = info["EOM"];
        int settlement_days = boost::lexical_cast<int>(settle.substr(0, settle.size() - 1));
            
        const Size size = 120;
        boost::shared_ptr<LiborForwardModelProcess> process(new LiborForwardModelProcess(size, index_));
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

                boost::shared_ptr<CalibrationHelper> swaptionHelper(new SwaptionHelper(expiry, tenor, swaptionVol, index_,
                                        index_->tenor(), parse_daycounter(fix_dct),
                                        index_->dayCounter(),
                                        getDiscountTermStructure(),
                                        CalibrationHelper::ImpliedVolError));

                swaptionHelper->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(LMMmodel_, getDiscountTermStructure())));
                modelCalibrator_.push_back(swaptionHelper); 
                CalibrationReport instr(key, instrument.expiry, instrument.tenor);
                calibrationReports_.push_back(instr);
            }
        }
        
    }
}

void LiborMarketModel::buildModel(){
    LevenbergMarquardt om(1e-6, 1e-6, 1e-6);
    LMMmodel_->calibrate(modelCalibrator_, om, EndCriteria(2000, 100, 1e-6, 1e-6, 1e-6));
    // measure the calibration error
    Real calculated = 0.0;
    for (int i=0; i<modelCalibrator_.size(); ++i) {
        calibrationReports_[i].diff = modelCalibrator_[i]->calibrationError();
        calibrationReports_[i].calculated = modelCalibrator_[i]->modelValue();
        calibrationReports_[i].expected = modelCalibrator_[i]->marketValue();
        calculated += std::pow(calibrationReports_[i].diff, 2.0);
    } 
}

void LiborMarketModel::generateMonteCarloPaths(int& numPaths, int& steps, double& maturity){
    const Size size = steps;
    Array parameters = LMMmodel_->params();
    Size length = parameters.size();
    boost::shared_ptr<LiborForwardModelProcess> process(new LiborForwardModelProcess(size, index_));

    boost::shared_ptr<LmCorrelationModel> corrModel(new LmLinearExponentialCorrelationModel(size, parameters[length - 2], parameters[length - 1]));

    boost::shared_ptr<LmVolatilityModel> volaModel(new LmLinearExponentialVolatilityModel(process->fixingTimes(), parameters[0], parameters[1], parameters[2], parameters[3]));

   // set-up pricing engine
    process->setCovarParam(boost::shared_ptr<LfmCovarianceParameterization>(new LfmCovarianceProxy(volaModel, corrModel)));

    // set-up a small Monte-Carlo simulation to price swations

    typedef PseudoRandom::rsg_type rsg_type;
    typedef MultiPathGenerator<rsg_type>::sample_type sample_type;

    std::vector<Time> tmp = process->fixingTimes();
    TimeGrid grid(tmp.begin(), tmp.end(), steps);
    timeGrid_.resize(steps);
    std::copy(grid.begin()+1, grid.end(), timeGrid_.begin());
    paths_.resize(numPaths*size/2, steps);

    Size i;
    std::vector<Size> location;
    for (i=0; i < tmp.size(); ++i) {
        location.push_back(std::find(grid.begin(),grid.end(),tmp[i])-grid.begin());
    }

    rsg_type rsg = PseudoRandom::make_sequence_generator(process->factors()*(grid.size()-1), BigNatural(42));

    const Size nrTrails = numPaths;
    MultiPathGenerator<rsg_type> generator(process, grid, rsg, false);

    boost::shared_ptr<LiborForwardModel>liborModel(new LiborForwardModel(process, volaModel, corrModel));

    Calendar calendar = index_->fixingCalendar();
    DayCounter dayCounter = index_->forwardingTermStructure()->dayCounter();
    BusinessDayConvention convention = index_->businessDayConvention();

    Date settlement  = index_->forwardingTermStructure()->referenceDate();
    int cnt = 0;
    boost::shared_ptr<SwaptionVolatilityMatrix> m = liborModel->getSwaptionVolatilityMatrix();
    for (i=1; i < size; ++i) {
        for (Size j=1; j <= size-i; ++j) {
            Date fwdStart    = settlement + Period(6*i, Months);
            Date fwdMaturity = fwdStart + Period(6*j, Months);


            if (i == j && i<=size/2) {

                GeneralStatistics stat;

                for (Size n=0; n<nrTrails; ++n) {
                    sample_type path = (n%2) ? generator.antithetic() : generator.next();

                    std::vector<Rate> rates(size);
                    for (Size k=0; k<process->size(); ++k) {
                        paths_(cnt, k) = path.value[k][location[i]];
                    }
                    cnt++;
                }
            }
        }
    }
}