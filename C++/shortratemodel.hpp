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
        double findb(double& dT, double& X);
        std::vector<double> solvedT(double& dT, double& X);
        void remove_duplicates(std::vector<double>& vec);
    public:
        CoxIngersollRossModel(std::string &modelDate, std::string &currency, std::string &exchange, std::string& interpolationType);
        CoxIngersollRossModel(std::string &modelDate, std::string &currency, std::string &exchange, std::string& interpolationType, double& theta, double& k, double& sigma, double& x0) : BaseModel(modelDate, currency, exchange, interpolationType) { CIRmodel_ = boost::shared_ptr<ExtendedCoxIngersollRoss>(new ExtendedCoxIngersollRoss(getDiscountTermStructure(), theta, k, sigma, x0)); 
            std::cout << CIRmodel_.get()->params() << std::endl; };
        void buildModel();
        std::vector<CalibrationReport> get_calibrationReport() { return calibrationReports_; };
        ProbabilityTree getLattice(std::vector<Period>& maturities);
};

class BlackKarasinskiModel : public BaseModel{
    private:
        std::set<std::string> fittable_;
        std::vector<boost::shared_ptr<CalibrationHelper>> modelCalibrator_;
        boost::shared_ptr<BlackKarasinski> BKmodel_;
        // For calibration;
        std::vector<CalibrationReport> calibrationReports_;
        // For lattice generation
        void buildProbabililtyTree(ProbabilityTree& ptree);
        double fsolve(double& init, double& targetPrice , boost::numeric::ublas::vector<double>& branches, boost::numeric::ublas::vector<double>& treeLevel, double& rateDelta, double& timeDelta);

        void parse_fitting_bucket(std::string& key, json& instrument);
        void parse_vol_instruments(std::string& key, json& instrument);
    public:
        BlackKarasinskiModel(std::string &modelDate, std::string &currency, std::string &exchange, std::string& interpolationType, double& a, double& sigma) : BaseModel(modelDate, currency, exchange, interpolationType) { BKmodel_ = boost::shared_ptr<BlackKarasinski>(new BlackKarasinski(getDiscountTermStructure(), a, sigma)); }
        BlackKarasinskiModel(std::string &modelDate, std::string &currency, std::string &exchange, std::string& interpolationType);
        void buildModel();
        ProbabilityTree getLattice(std::vector<Period>& maturities);
        std::vector<CalibrationReport> getCalibrationReport() { return calibrationReports_; };
};

#endif