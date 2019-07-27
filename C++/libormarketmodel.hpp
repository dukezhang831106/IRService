#ifndef libor_market_model_hpp
#define libor_market_model_hpp

#include <ql/time/date.hpp>
#include <ql/time/calendars/all.hpp>
#include <ql/time/calendars/jointcalendar.hpp>
#include <ql/time/daycounters/actual360.hpp>

#include <ql/indexes/ibor/euribor.hpp>

#include <ql/termstructures/yield/zerocurve.hpp>

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
        //std::vector<> yieldCurve_;
        void parse_calendar(std::string& holmasks);
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