#ifndef libor_market_model_hpp
#define libor_market_model_hpp

#include <ql/time/date.hpp>
#include <ql/time/calendars/all.hpp>
#include <ql/time/calendars/jointcalendar.hpp>
#include <ql/time/daycounters/actual360.hpp>

#include <ql/indexes/ibor/euribor.hpp>
#include <ql/indexes/ibor/usdlibor.hpp>

#include <ql/termstructures/yield/zerocurve.hpp>
#include <ql/termstructures/yield/ratehelpers.hpp>

#include <ql/models/shortrate/calibrationhelpers/caphelper.hpp>
#include <ql/models/shortrate/calibrationhelpers/swaptionhelper.hpp>

#include <ql/legacy/libormarketmodels/lmlinexpcorrmodel.hpp>
#include <ql/legacy/libormarketmodels/lmextlinexpvolmodel.hpp>
#include <ql/legacy/libormarketmodels/liborforwardmodel.hpp>
#include <ql/legacy/libormarketmodels/lfmswaptionengine.hpp>

#include <ql/math/optimization/levenbergmarquardt.hpp>

#include <ql/pricingengines/capfloor/analyticcapfloorengine.hpp>

#include <ql/quotes/simplequote.hpp>
#include <ql/utilities/dataparsers.hpp>
#include <ql/time/businessdayconvention.hpp>
#include <ql/time/daycounter.hpp>
#include <ql/time/period.hpp>



#include <pqxx/pqxx>

#include <boost/lexical_cast.hpp>

#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include "json.hpp"

using namespace QuantLib;
using json = nlohmann::json;

struct GeneralInstrumentInformation{
    Date date;
    std::string currency, instrumentType, expiry, tenor, source;
    double quote;
    GeneralInstrumentInformation(pqxx::row& record) {        
        for(int i = 0; i < record.size(); i++){
            pqxx::field f(record, i);
            if(i == 0){
                Date d = DateParser::parseFormatted(f.c_str(), "%Y-%m-%d");
            }
            else if (i == 1){
                currency = f.c_str();
            }
            else if (i == 2){
                instrumentType = f.c_str();
            }
            else if (i == 3){
                expiry = f.c_str();
            }
            else if (i == 4){
                tenor = f.c_str();
            }
            else if (i == 5){
                source = f.c_str();
            }
            else if (i == 6){
                quote = boost::lexical_cast<double>(f.c_str());
            }
        }
    }
};


class LiborMarketModel {
    private:
        std::string currency_, exchange_, modelModelDate_;
        Calendar calendar_;
        std::vector<boost::shared_ptr<CalibrationHelper>> curveCalibrator_;
        std::vector<boost::shared_ptr<CalibrationHelper>> modelCalibrator_;
        RelinkableHandle<YieldTermStructure> termStructure_;
    
        Calendar parse_calendar(std::string& holmasks);
        BusinessDayConvention parse_buisnessdayconvention(std::string& bdc);
        DayCounter parse_daycounter(std::string& dct);
        Period parse_period(std::string& p);
        void parse_instruments(std::string& key, json& instrument);
        ext::shared_ptr<IborIndex> makeIndex(std::vector<Date> dates, std::vector<Rate> rates);
        ext::shared_ptr<IborIndex> makeIndex();

    public:
        LiborMarketModel(std::string &modelDate, std::string &currency, std::string &exchange);
        void setCurveCalibrator();
        void setModelCalibrator();
        void buildCurve();
        void buildModel();
        void test_connection();
};

#endif