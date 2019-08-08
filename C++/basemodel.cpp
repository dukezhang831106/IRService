#include "basemodel.hpp"

BaseModel::BaseModel(std::string& modelDate, std::string& currency, std::string& exchange, std::string& interpolationType) : modelDate_(modelDate), modelModelDate_(DateParser::parseISO(modelDate)), currency_(currency), exchange_(exchange), interpolationType_(interpolationType) {
    std::ifstream file("C++/data/configuration.json");
    json output;
    file >> output;
    for(auto & exchange : output[currency_]){
        if (exchange["Exchange"] == exchange_){
            // parse general information
            std::string settle = exchange["SettlementDays"], daycounter = exchange["DayCountConvention"], holmask = exchange["Holidays"], dayCounter_ = exchange["DayCountConvention"];
            calendar_ = parse_calendar(holmask);
            int fixingDays = boost::lexical_cast<int>(settle.substr(0, settle.size() - 1));

            // set up model date
            Date settlementDate = calendar_.advance(modelModelDate_, fixingDays, Days);
            // settlement must be a business day
            settlementDate_ = calendar_.adjust(settlementDate);
            modelModelDate_ = calendar_.advance(settlementDate, -fixingDays, Days);
            // nothing to do with Date::todaysDate
            Settings::instance().evaluationDate() = modelModelDate_;
            modelModelDate_ = Settings::instance().evaluationDate();

            for (auto & instrument : exchange["CurveInstruments"].items()){
                std::string key(instrument.key());
                parse_linear_instruments(key, instrument.value());
            }
            buildCurve();
        }
    }
    return;
}

void BaseModel::parse_linear_instruments(std::string& key, json& info){
    pqxx::connection conn("dbname=EikonInstrument user=postgres host=127.0.0.1 port=5432 password=123456");
    pqxx::work txn(conn);

    std::stringstream command;
    command << "SELECT * " << "FROM \"IRInstruments\" " << "WHERE \"QuoteDate\"=Date('" 
        << modelDate_ << "') AND \"Currency\"='USD' AND \"InstrumentType\"='" << key << "';";
    pqxx::result res = txn.exec(command.str());
    
    if (key == "Depo"){
        std::string dct = info["DayCountConvention"], bdc = info["BusinessDayConvention"], settle = info["SettlementDays"];
        bool eom = info["EOM"];
        int settlement_days = boost::lexical_cast<int>(settle.substr(0, settle.size() - 1));
        for(int i = 0; i < res.size(); i++){
            pqxx::row record = res[i];
            GeneralInstrumentInformation instrument(record);
            boost::shared_ptr<Quote> rate(new SimpleQuote(instrument.quote/100.0));
            boost::shared_ptr<RateHelper> holder(new DepositRateHelper(Handle<Quote>(rate), parse_period(instrument.expiry), settlement_days, calendar_, parse_buisnessdayconvention(bdc), eom, parse_daycounter(dct)));
            curveData_.insert(std::make_pair(parse_period(instrument.expiry), instrument.quote/100.0));
            curveCalibrator_.push_back(holder);
        }
    }
    else if (key == "Fut"){
        std::string dct = info["DayCountConvention"], bdc = info["BusinessDayConvention"], settle = info["SettlementDays"], freq = info["Frequency"];
        bool eom = info["EOM"];
        int settlement_days = boost::lexical_cast<int>(settle.substr(0, settle.size() - 1)), frequency = boost::lexical_cast<int>(freq.substr(0, freq.size() - 1));
        for(int i = 0; i < res.size(); i++){
            pqxx::row record = res[i];
            GeneralInstrumentInformation instrument(record);
            boost::shared_ptr<Quote> price(new SimpleQuote(instrument.quote));
            boost::shared_ptr<RateHelper> holder(new FuturesRateHelper(Handle<Quote>(price), IMM::date(instrument.expiry), frequency, calendar_, parse_buisnessdayconvention(bdc), eom, parse_daycounter(dct)));
            curveCalibrator_.push_back(holder);
        }
    }
    else if (key == "Swap"){
        std::string fix_dct = info["FixedLegDayCountConvention"], fix_bdc = info["FixedLegBusinessDayConvention"], fix_freq = info["FixedLegFrequency"];
        std::string flt_dct = info["FloatingLegDayCountConvention"], flt_bdc = info["FloatingLegBusinessDayConvention"], flt_freq = info["FloatingLegFrequency"];
        std::string settle = info["SettlementDays"], indices = info["Index"];
        indices_ = indices;
        bool eom = info["EOM"];
        int settlement_days = boost::lexical_cast<int>(settle.substr(0, settle.size() - 1));
        for(int i = 0; i < res.size(); i++){
            pqxx::row record = res[i];
            GeneralInstrumentInformation instrument(record);
            boost::shared_ptr<Quote> rate(new SimpleQuote(instrument.quote/100.0));
            boost::shared_ptr<RateHelper> holder(new SwapRateHelper(Handle<Quote>(rate), parse_period(instrument.expiry), calendar_, 
                parse_frequency(fix_freq), parse_buisnessdayconvention(fix_bdc), parse_daycounter(fix_dct), parse_indices(indices)));
            curveData_.insert(std::make_pair(parse_period(instrument.expiry), instrument.quote/100.0));
            curveCalibrator_.push_back(holder);
        }
    }
}

BusinessDayConvention BaseModel::parse_buisnessdayconvention(std::string& bdc){
    if (bdc == "ModifiedFollwing")
        return ModifiedFollowing;
    else if (bdc == "Unadjusted")
        return Unadjusted;
    return Unadjusted;
}

DayCounter BaseModel::parse_daycounter(std::string& dct){
    if (dct == "Actual/360")
        return Actual360();
    else if (dct == "30/360")
        return Thirty360();
    else if (dct == "Actual/365F")
        return Actual365Fixed();
    else if (dct == "Actual/Actual")
        return ActualActual();
    return Actual360();
}

Frequency BaseModel::parse_frequency(std::string& freq){
    if (freq == "1M")
        return Monthly;
    else if (freq == "3M")
        return Quarterly;
    else if (freq == "6M")
        return Semiannual;
    else if (freq == "1Y" || freq == "12M")
        return Annual;
    return Semiannual;
}

boost::shared_ptr<IborIndex> BaseModel::parse_indices(std::string& indices){
    std::string ccy(indices.substr(0, 3)), tenor(indices.substr(3, indices.size()));
    if (ccy == "USD"){
        if (tenor == "3M"){
            return boost::shared_ptr<IborIndex>(new USDLibor(parse_period(tenor)));
        }
    }
}

Period BaseModel::parse_period(std::string& p){
    if (p == "TN")
        return Period(1, Days);
    else if (p == "SN")
        return Period(2, Days);
    else if (p[p.size()-1] == 'D'){
        int length = boost::lexical_cast<int>(p.substr(0, p.size()-1));
        return Period(length, Days);
    }
    else if (p[p.size()-1] == 'M'){
        int length = boost::lexical_cast<int>(p.substr(0, p.size()-1));
        return Period(length, Months);
    }
        
    else if (p[p.size()-1] == 'Y'){
        int length = boost::lexical_cast<int>(p.substr(0, p.size()-1));
        return Period(length, Years);
    }
    return Period(0, Days);
}

Calendar BaseModel::parse_calendar(std::string& holmasks) {
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

void BaseModel::buildCurve(){
    double tolerance = 1e-10;
    // First we build a discount curve
    auto oisIndex = boost::make_shared<FedFunds>();

    // create container for rate helpers
    std::vector<boost::shared_ptr<RateHelper>> rateHelpers;

    Natural settlementDays = settlementDate_ - modelModelDate_;
    // create first 1d cash instrument for eonia curve - using deposit rate helper
    auto Q1D = boost::make_shared<SimpleQuote>(curveCalibrator_[0].get()->quote().currentLink().get()->value());
    rateHelpers.push_back(boost::make_shared<DepositRateHelper>(Handle<Quote>(Q1D), Period(1, Days), oisIndex->fixingDays(), oisIndex->fixingCalendar(), oisIndex->businessDayConvention(), oisIndex->endOfMonth(), oisIndex->dayCounter()));

    // create source data for eonia curve (period, rate)
    

    // create other instruments for eonia curve - using ois rate helper
    std::for_each(curveData_.begin(), curveData_.end(), [settlementDays, &rateHelpers, &oisIndex](std::pair<Period, Real> p) -> void  { rateHelpers.push_back(boost::make_shared<OISRateHelper>(settlementDays, p.first, Handle<Quote>(boost::make_shared<SimpleQuote>(p.second)), oisIndex)); });

    // create piecewise term structure
    if (interpolationType_ == "LogLinear"){
        boost::shared_ptr<YieldTermStructure> oisCurve(boost::make_shared<PiecewiseYieldCurve<Discount, LogLinear>>(settlementDays, oisIndex->fixingCalendar(), rateHelpers, oisIndex->dayCounter()));
        discountingTermStructure_.linkTo(oisCurve);
        boost::shared_ptr<YieldTermStructure> fwdCurve(new PiecewiseYieldCurve<ZeroYield, LogLinear>(settlementDate_, curveCalibrator_, parse_daycounter(dayCounter_), tolerance));
        forecastingTermStructure_.linkTo(fwdCurve);
        return;
    }
    else if(interpolationType_ == "Linear"){
        boost::shared_ptr<YieldTermStructure> oisCurve(boost::make_shared<PiecewiseYieldCurve<Discount, Linear>>(settlementDays, oisIndex->fixingCalendar(), rateHelpers, oisIndex->dayCounter()));
        discountingTermStructure_.linkTo(oisCurve);
        boost::shared_ptr<YieldTermStructure> fwdCurve(new PiecewiseYieldCurve<ZeroYield, Linear>(settlementDate_, curveCalibrator_, parse_daycounter(dayCounter_), tolerance));
        forecastingTermStructure_.linkTo(fwdCurve);
        return;
    }
    else if(interpolationType_ == "Cubic"){
        boost::shared_ptr<YieldTermStructure> oisCurve(boost::make_shared<PiecewiseYieldCurve<Discount, Cubic>>(settlementDays, oisIndex->fixingCalendar(), rateHelpers, oisIndex->dayCounter()));
        discountingTermStructure_.linkTo(oisCurve);
        boost::shared_ptr<YieldTermStructure> fwdCurve(new PiecewiseYieldCurve<ZeroYield, Cubic>(settlementDate_, curveCalibrator_, parse_daycounter(dayCounter_), tolerance));
        forecastingTermStructure_.linkTo(fwdCurve);
        return;
    }
}

boost::shared_ptr<IborIndex> BaseModel::makeIndex(std::string& indices) {
    std::string ccy(indices.substr(0, 3)), tenor(indices.substr(3, indices.size()));
    if (ccy == "USD"){
        if (tenor == "3M"){
            RelinkableHandle<YieldTermStructure> termStructure = getDiscountTermStructure();
            boost::shared_ptr<IborIndex> index(new USDLibor(parse_period(tenor), termStructure));
            
            return index;
        }
    }
    return boost::shared_ptr<IborIndex>(new USDLibor(parse_period(tenor), getDiscountTermStructure()));
}

void BaseModel::getCurveZeroRates(std::vector<Time>& times){    
    for(int i = 0; i < times.size(); i++){
        Rate r = discountingTermStructure_.currentLink().get()->zeroRate(times[i], Compounded, Semiannual);
        std::cout << times[i] << ": " << r << std::endl;
    }
}

void BaseModel::getCashDF(std::vector<Time>& times){    
    for(int i = 0; i < times.size(); i++){
        Rate r = discountingTermStructure_.currentLink().get()->discount(times[i]);
        std::cout << times[i] << ": " << r << std::endl;
    }
}