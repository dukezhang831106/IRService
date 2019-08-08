#ifndef short_rate_model_hpp
#define short_rate_model_hpp

#include "basemodel.hpp"

class HullWhiteModel : public BaseModel{
    private:
        std::set<std::string> fittable_;
        std::vector<boost::shared_ptr<CalibrationHelper>> modelCalibrator_;
        boost::shared_ptr<HullWhite> HWmodel_;
        // For calibration;
        std::vector<CalibrationReport> calibrationReports_;
        // For Monte Carlo Simulation
        boost::numeric::ublas::vector<double> timeGrid_;
        boost::numeric::ublas::matrix<double> paths_;
        void parse_fitting_bucket(std::string& key, json& instrument);
        void parse_vol_instruments(std::string& key, json& instrument);

    public:
        HullWhiteModel(std::string &modelDate, std::string &currency, std::string &exchange, std::string& interpolationType);
        void buildModel();
        std::vector<CalibrationReport> get_calibrationReport() { return calibrationReports_; };
        void generateMonteCarloPaths(int& numPaths, int& steps, double& maturity);
        boost::numeric::ublas::vector<double> getTimeGrid() { return timeGrid_; };
        boost::numeric::ublas::matrix<double> getMonteCarloPaths() { return paths_; };
};


class CoxIngersollRossModel : public BaseModel{
    private:
        std::set<std::string> fittable_;
        std::vector<boost::shared_ptr<CalibrationHelper>> modelCalibrator_;
        boost::shared_ptr<ExtendedCoxIngersollRoss> CIRmodel_;
        // For calibration;
        std::vector<CalibrationReport> calibrationReports_;
        // For Monte Carlo Simulation
        boost::numeric::ublas::vector<double> timeGrid_;
        boost::numeric::ublas::matrix<double> paths_;
        void parse_fitting_bucket(std::string& key, json& instrument);
        void parse_vol_instruments(std::string& key, json& instrument);

    public:
        CoxIngersollRossModel(std::string &modelDate, std::string &currency, std::string &exchange, std::string& interpolationType);
        void buildModel();
        std::vector<CalibrationReport> get_calibrationReport() { return calibrationReports_; };
        void generateMonteCarloPaths(int& nodes, int& steps, double& maturity);
        boost::numeric::ublas::vector<double> getTimeGrid() { return timeGrid_; };
        boost::numeric::ublas::matrix<double> getMonteCarloPaths() { return paths_; };
};

class BlackKarasinskiModel : public BaseModel{
    private:
        std::set<std::string> fittable_;
        std::vector<boost::shared_ptr<CalibrationHelper>> modelCalibrator_;
        boost::shared_ptr<BlackKarasinski> BKmodel_;
        // For calibration;
        std::vector<CalibrationReport> calibrationReports_;
        // For Monte Carlo Simulation
        boost::numeric::ublas::vector<double> timeGrid_;
        boost::numeric::ublas::matrix<double> paths_;
        void parse_fitting_bucket(std::string& key, json& instrument);
        void parse_vol_instruments(std::string& key, json& instrument);
    public:
        BlackKarasinskiModel(std::string &modelDate, std::string &currency, std::string &exchange, std::string& interpolationType);
        void buildModel();
        std::vector<CalibrationReport> get_calibrationReport() { return calibrationReports_; };
        void generateMonteCarloPaths(int& numPaths, int& steps, double& maturity);
        boost::numeric::ublas::vector<double> getTimeGrid() { return timeGrid_; };
        boost::numeric::ublas::matrix<double> getMonteCarloPaths() { return paths_; };
};


#endif