/* General model holder */
/* Author: Tianyi Zhang */
/* support Interest Rate Curve and Interest Rate Model */
#ifndef base_model_hpp
#define base_model_hpp

 // avoid pragma warning
#define BOOST_PENDING_INTEGER_LOG2_HPP
#include <boost/integer/integer_log2.hpp>

#include <ql/quantlib.hpp>

#include <pqxx/pqxx>

#include <boost/lexical_cast.hpp>

#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <set>
#include "json.hpp"

using namespace QuantLib;
using json = nlohmann::json;

struct GeneralInstrumentInformation{
    Date date;
    std::string currency, instrumentType, expiry, tenor, source;
    double quote;
    GeneralInstrumentInformation(pqxx::row& record) {        
        for(int i = 0; i < record.size(); i++){
            pqxx::field f = record[i];
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
                source = f.c_str();
            }
            else if (i == 4){
                expiry = f.c_str();
            }
            else if (i == 5){
                tenor = f.c_str();
            }
            else if (i == 6){
                quote = boost::lexical_cast<double>(f.c_str());
            }
        }
    }
};

struct CalibrationReport{
    std::string instrumentType, expiry, tenor, source;
    double expected, calculated, diff;
};

class BaseModel {
    private:
        std::string currency_, exchange_, modelDate_, interpolationType_, indices_, dayCounter_;
        Date modelModelDate_, settlementDate_;
        Calendar calendar_;
        std::map<Period, Real> curveData_;
        std::vector<boost::shared_ptr<RateHelper>> curveCalibrator_;
        RelinkableHandle<YieldTermStructure> discountingTermStructure_, forecastingTermStructure_;
    public:
        boost::shared_ptr<IborIndex> makeIndex(std::string& indices);
        BaseModel(std::string &modelDate, std::string &currency, std::string &exchange, std::string& interpolationType);
        Calendar parse_calendar(std::string& holmasks);
        BusinessDayConvention parse_buisnessdayconvention(std::string& bdc);
        DayCounter parse_daycounter(std::string& dct);
        Period parse_period(std::string& p);
        Frequency parse_frequency(std::string& freq);
        void parse_linear_instruments(std::string& key, json& instrument);

        boost::shared_ptr<IborIndex> parse_indices(std::string& indices);
        void buildCurve();
        virtual void buildModel() = 0;
        std::vector<CalibrationReport> generateCurveCalibrationReport();
        //virtual std::vector<CalibrationReport> generateModelCalibrationReport() = 0;
        void getCurveZeroRates(std::vector<Time>& times);
        void getCashDF(std::vector<Time>& times);

        std::string getCurrency() { return currency_; }
        std::string getExchange() { return exchange_; }
        std::string getModelDate() { return modelDate_; }
        Date getModelModelDate() { return modelModelDate_;}
        Date getSettlementDate() { return settlementDate_; }
        Calendar getCalendar() { return calendar_; }
        std::vector<boost::shared_ptr<RateHelper>> getCurveCalibrator() { return curveCalibrator_; }
        RelinkableHandle<YieldTermStructure> getDiscountTermStructure() { return discountingTermStructure_; }
        RelinkableHandle<YieldTermStructure> getForecastTermStructure() { return forecastingTermStructure_; }

};

#endif