#include "libormarketmodel.hpp"

LiborMarketModel::LiborMarketModel(std::string& modelDate, std::string& currency, std::string& exchange) : modelModelDate_(modelDate), currency_(currency), exchange_(exchange) {
    std::ifstream file("C++/data/configuration.json");
    json output;
    file >> output;
    for(auto & exchange : output[currency_]){
        if (exchange["Exchange"] == exchange_){
            std::string holmask = exchange["Holidays"];
            calendar_ = parse_calendar(holmask);
            
            for (auto & instrument : exchange["CurveInstruments"].items()){
                std::string key(instrument.key());
                parse_instruments(key, instrument.value());
            }

            for (auto & vol_instrument : exchange["VolatilityInstruments"].items()){
                std::string key(vol_instrument.key());
                parse_instruments(key, vol_instrument.value());
            }
        }
    }
    return;    
}

void LiborMarketModel::parse_instruments(std::string& key, json& info){
    pqxx::connection conn("dbname=EikonInstrument user=postgres host=127.0.0.1 port=1111 password=123456");
    pqxx::work txn(conn);

    std::stringstream command;
    command << "SELECT * " << "FROM \"IRInstruments\" " << "WHERE \"QuoteDate\"=Date('" 
        << modelModelDate_ << "') AND \"Currency\"='USD' AND \"InstrumentType\"='" << key << "';";
    pqxx::result res = txn.exec(command.str());
    
    if (key == "Depo"){
        std::string dct = info["DayCountConvention"], bdc = info["BusinessDayConvention"];
        
        for(int i = 0; i < res.size(); i++){
            pqxx::row record(res, i);
            GeneralInstrumentInformation instrument(record);
            Handle<Quote> r(instrument.quote);
            curveCalibrator_.push_back(boost::shared_ptr<CalibrationHelper>(new DepositRateHelper(r, parse_period(instrument.expiry), 
                2, calendar_, bdc, false, dct)));
            std::cout << instrument.quote << std::endl;
            /*ext::shared_ptr<BlackCalibrationHelper> depohelper(
            new DepositRateHelper(, capVol, index, Annual,
                          index->dayCounter(), true, termStructure,
                          BlackCalibrationHelper::ImpliedVolError));*/
        }
        
    }
}

BusinessDayConvention LiborMarketModel::parse_buisnessdayconvention(std::string& bdc){
    if (bdc == "ModifiedFollwing")
        return ModifiedFollowing;
    return Unadjusted;
}


DayCounter LiborMarketModel::parse_daycounter(std::string& dct){
    if (dct == "Actual/360")
        return Actual360();
    return Actual360();
}

Period LiborMarketModel::parse_period(std::string& p){
    int length = boost::lexical_cast<int>(p.substr(0, p.size()-1));
    if (p == "TN")
        return Period(1, Days);
    else if (p == "SN")
        return Period(2, Days);
    else if (p[p.size()-1] == 'D')
        return Period(length, Days);
    else if (p[p.size()-1] == 'M')
        return Period(length, Months);
    else if (p[p.size()-1] == 'Y')
        return Period(length, Years);
    return Period(length, Years);
}

Calendar LiborMarketModel::parse_calendar(std::string& holmasks) {
    std::stringstream holidays(holmasks);
    int cnt = 0;
    Calendar calendar;
    while(holidays.good()){
        std::string hol;
        getline(holidays, hol, ',');
        Calendar cal;
        if (hol == "NY") 
            cal = UnitedStates();
        else if (hol == "LON")
            cal = UnitedKingdom();
        else if (hol == "EUR")
            cal = TARGET();
        else if (hol == "TOK")
            cal = Japan();
        else if (hol == "BEI")
            cal = China();
        else if (hol == "SYD")
            cal = Australia();
        else if (hol == "ZUR")
            cal = Switzerland();
        else 
            cal = UnitedStates();
        if (cnt == 0){
            calendar = cal;
        }
        else{
            calendar = JointCalendar(calendar, cal);
        }
        cnt++;
    }
    return calendar;
}

void LiborMarketModel::buildCurve(){
    std::ifstream file("C++/data/configuration.json");
    json output;
    file >> output;
    for(auto & exchange : output[currency_]){
        if (exchange["Exchange"] == exchange_){
            //Calendar cal = parse_calendar(std::string(exchange["Holiday"]));
            std::cout << exchange["Instruments"] << std::endl;
        }
    }
    return;
}

ext::shared_ptr<IborIndex> LiborMarketModel::makeIndex(std::vector<Date> dates, std::vector<Rate> rates) {
    DayCounter dayCounter = Actual360();
    RelinkableHandle<YieldTermStructure> termStructure;
    ext::shared_ptr<IborIndex> index(new Euribor6M(termStructure));
    Date todaysDate = index->fixingCalendar().adjust(Date(4,September,2005));
    Settings::instance().evaluationDate() = todaysDate;

    dates[0] = index->fixingCalendar().advance(todaysDate, index->fixingDays(), Days);

    termStructure.linkTo(ext::shared_ptr<YieldTermStructure>(new ZeroCurve(dates, rates, dayCounter)));

    return index;
}

ext::shared_ptr<IborIndex> LiborMarketModel::makeIndex() {
    std::vector<Date> dates;
    std::vector<Rate> rates;
    dates.push_back(Date(4,September,2005));
    dates.push_back(Date(4,September,2018));
    rates.push_back(0.039);
    rates.push_back(0.041);

    return makeIndex(dates, rates);
}

void LiborMarketModel::buildModel(){
    const Size size = 14;
    const Real tolerance = 1e-3;

    Volatility capVols[] = {0.145708,0.158465,0.166248,0.168672,
                            0.169007,0.167956,0.166261,0.164239,
                            0.162082,0.159923,0.157781,0.155745,
                            0.153776,0.151950,0.150189,0.148582,
                            0.147034,0.145598,0.144248};

    Volatility swaptionVols[] = {0.170595, 0.166844, 0.158306, 0.147444,
                                 0.136930, 0.126833, 0.118135, 0.175963,
                                 0.166359, 0.155203, 0.143712, 0.132769,
                                 0.122947, 0.114310, 0.174455, 0.162265,
                                 0.150539, 0.138734, 0.128215, 0.118470,
                                 0.110540, 0.169780, 0.156860, 0.144821,
                                 0.133537, 0.123167, 0.114363, 0.106500,
                                 0.164521, 0.151223, 0.139670, 0.128632,
                                 0.119123, 0.110330, 0.103114, 0.158956,
                                 0.146036, 0.134555, 0.124393, 0.115038,
                                 0.106996, 0.100064};

    ext::shared_ptr<IborIndex> index = makeIndex();
    ext::shared_ptr<LiborForwardModelProcess> process(
        new LiborForwardModelProcess(size, index));
    Handle<YieldTermStructure> termStructure = index->forwardingTermStructure();

    // set-up the model
    ext::shared_ptr<LmVolatilityModel> volaModel(
                    new LmExtLinearExponentialVolModel(process->fixingTimes(),
                                                       0.5,0.6,0.1,0.1));

    ext::shared_ptr<LmCorrelationModel> corrModel(
                     new LmLinearExponentialCorrelationModel(size, 0.5, 0.8));

    ext::shared_ptr<LiborForwardModel> model(
                        new LiborForwardModel(process, volaModel, corrModel));

    Size swapVolIndex = 0;
    DayCounter dayCounter=index->forwardingTermStructure()->dayCounter();

    // set-up calibration helper
    std::vector<ext::shared_ptr<BlackCalibrationHelper> > calibrationHelper;

    Size i;
    for (i=2; i < size; ++i) {
        const Period maturity = i*index->tenor();
        Handle<Quote> capVol(
            ext::shared_ptr<Quote>(new SimpleQuote(capVols[i-2])));

        ext::shared_ptr<BlackCalibrationHelper> caphelper(
            new CapHelper(maturity, capVol, index, Annual,
                          index->dayCounter(), true, termStructure,
                          BlackCalibrationHelper::ImpliedVolError));

        caphelper->setPricingEngine(ext::shared_ptr<PricingEngine>(
                           new AnalyticCapFloorEngine(model, termStructure)));

        calibrationHelper.push_back(caphelper);

        if (i<= size/2) {
            // add a few swaptions to test swaption calibration as well
            for (Size j=1; j <= size/2; ++j) {
                const Period len = j*index->tenor();
                Handle<Quote> swaptionVol(
                    ext::shared_ptr<Quote>(
                        new SimpleQuote(swaptionVols[swapVolIndex++])));

                ext::shared_ptr<BlackCalibrationHelper> swaptionHelper(
                    new SwaptionHelper(maturity, len, swaptionVol, index,
                                       index->tenor(), dayCounter,
                                       index->dayCounter(),
                                       termStructure,
                                       BlackCalibrationHelper::ImpliedVolError));

                swaptionHelper->setPricingEngine(
                     ext::shared_ptr<PricingEngine>(
                                 new LfmSwaptionEngine(model,termStructure)));

                calibrationHelper.push_back(swaptionHelper);
            }
        }
    }

    LevenbergMarquardt om(1e-6, 1e-6, 1e-6);
    model->calibrate(calibrationHelper, om, EndCriteria(2000, 100, 1e-6, 1e-6, 1e-6));

    // measure the calibration error
    Real calculated = 0.0;
    for (i=0; i<calibrationHelper.size(); ++i) {
        Real diff = calibrationHelper[i]->calibrationError();
        calculated += diff*diff;
    }

    if (std::sqrt(calculated) > tolerance){
        std::cout << "Failed to calibrate libor forward model" << std::endl;
        std::cout << "Calculated diff: " << std::sqrt(calculated) << std::endl;
        std::cout << "Expected : smaller than  " << tolerance << std::endl;
    }

    std::cout << "Success" << std::endl;
    std::cout << "Calculated diff: " << std::sqrt(calculated) << std::endl;
    std::cout << "Expected : smaller than  " << tolerance << std::endl;
}

/*void LiborMarketModel::test_connection(){
    pqxx::connection conn("dbname=EikonInstrument user=postgres host=127.0.0.1 port=1111 password=123456");
    pqxx::work txn(conn);

    std::stringstream command;
    command << "SELECT * " <<
        "FROM \"IRInstruments\" " <<
        "WHERE \"QuoteDate\"=Date('" << modelModelDate_ << "') AND \"Currency\"='USD' AND \"InstrumentType\"='Swaption';";

    pqxx::result res = txn.exec(command.str());
    
    for(int i = 0; i < res.size(); i++){
        for(int j = 0; j < res[i].size(); j++){
            std::cout << res[i][j] << ", ";
        }
        std::cout << std::endl;
    }
    return;
}*/