#ifndef libor_market_model_hpp
#define libor_market_model_hpp

 // avoid pragma warning
#include "basemodel.hpp"

class LiborMarketModel : public BaseModel{
    private:
        std::set<std::string> fittable_;
        std::vector<boost::shared_ptr<CalibrationHelper>> modelCalibrator_;
        boost::shared_ptr<LiborForwardModel> LMMmodel_;
        void parse_fitting_bucket(std::string& key, json& instrument);
        void parse_vol_instruments(std::string& key, json& instrument);

    public:
        LiborMarketModel(std::string &modelDate, std::string &currency, std::string &exchange, std::string& interpolationType) : BaseModel(modelDate, currency, exchange, interpolationType) {};
        void buildModel();
        std::vector<CalibrationReport> generateModelCalibrationReport();
};

#endif