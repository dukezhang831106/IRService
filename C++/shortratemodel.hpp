#ifndef short_rate_model_hpp
#define short_rate_model_hpp

 // avoid pragma warning
#include "basemodel.hpp"

class HullWhiteModel : public BaseModel{
    private:
        std::set<std::string> fittable_;
        std::vector<boost::shared_ptr<CalibrationHelper>> modelCalibrator_;
        boost::shared_ptr<HullWhite> HWmodel_;
        void parse_fitting_bucket(std::string& key, json& instrument);
        void parse_vol_instruments(std::string& key, json& instrument);

    public:
        HullWhiteModel(std::string &modelDate, std::string &currency, std::string &exchange, std::string& interpolationType);
        void buildModel();
        std::vector<CalibrationReport> generateModelCalibrationReport();
};


class ExtendedCIRModel : public BaseModel{
    private:
        std::set<std::string> fittable_;
        std::vector<boost::shared_ptr<CalibrationHelper>> modelCalibrator_;
        boost::shared_ptr<ExtendedCoxIngersollRoss> eCIRmodel_;
        void parse_fitting_bucket(std::string& key, json& instrument);
        void parse_vol_instruments(std::string& key, json& instrument);

    public:
        ExtendedCIRModel(std::string &modelDate, std::string &currency, std::string &exchange, std::string& interpolationType);
        void buildModel();
        std::vector<CalibrationReport> generateModelCalibrationReport();
};

class BlackKarasinskiModel : public BaseModel{
    private:
        std::set<std::string> fittable_;
        std::vector<boost::shared_ptr<CalibrationHelper>> modelCalibrator_;
        boost::shared_ptr<BlackKarasinski> BKmodel_;
        void parse_fitting_bucket(std::string& key, json& instrument);
        void parse_vol_instruments(std::string& key, json& instrument);

    public:
        BlackKarasinskiModel(std::string &modelDate, std::string &currency, std::string &exchange, std::string& interpolationType);
        void buildModel();
        std::vector<CalibrationReport> generateModelCalibrationReport();
};


#endif