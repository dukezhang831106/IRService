#ifndef libor_market_model_hpp
#define libor_market_model_hpp

#include <ql/time/date.hpp>
#include <ql/time/calendars/all.hpp>
#include <ql/time/calendars/jointcalendar.hpp>
#include <ql/models/calibrationhelper.hpp>

#include <ql/utilities/dataparsers.hpp>

#include <pqxx/pqxx>

#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include "json.hpp"

using namespace QuantLib;
using json = nlohmann::json;

class LiborMarketModel {
    private:
        std::string currency_, exchange_, modelModelDate_;
        Calendar calendar_;
        std::vector<boost::shared_ptr<CalibrationHelper>> curveCalibrator_;
        std::vector<boost::shared_ptr<CalibrationHelper>> modelCalibrator_;
        std::vector<> yieldCurve_;
        void parse_calendar(std::string& holmasks);

    public:
        LiborMarketModel(std::string modelDate, std::string currency = "USD", std::string exchange = "NY");
        void setCurveCalibrator();
        void setModelCalibrator();
        void buildCurve();
        void buildModel();
        void test_connection();
};

#endif