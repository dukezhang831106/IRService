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
        CoxIngersollRossModel(std::string &modelDate, std::string &currency, std::string &exchange, std::string& interpolationType, double& theta, double& k, double& sigma, double& x0) : BaseModel(modelDate, currency, exchange, interpolationType) { CIRmodel_ = boost::shared_ptr<ExtendedCoxIngersollRoss>(new ExtendedCoxIngersollRoss(getDiscountTermStructure(), theta, k, sigma, x0)); };
        void buildModel();
        std::vector<CalibrationReport> get_calibrationReport() { return calibrationReports_; };
        void generateMonteCarloPaths(int& nodes, int& steps, double& maturity);
        boost::numeric::ublas::vector<double> getTimeGrid() { return timeGrid_; };
        boost::numeric::ublas::matrix<double> getMonteCarloPaths() { return paths_; };
};

class BlackKarasinskiModel : public BaseModel{
    private:
        class BlackKarasinskiProcess;
        class BlackKarasinskiHelper;
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
        BlackKarasinskiModel(std::string &modelDate, std::string &currency, std::string &exchange, std::string& interpolationType, double& a, double& sigma) : BaseModel(modelDate, currency, exchange, interpolationType) { BKmodel_ = boost::shared_ptr<BlackKarasinski>(new BlackKarasinski(getDiscountTermStructure(), a, sigma)); }
        BlackKarasinskiModel(std::string &modelDate, std::string &currency, std::string &exchange, std::string& interpolationType);
        void buildModel();
        std::vector<CalibrationReport> get_calibrationReport() { return calibrationReports_; };
        void generateMonteCarloPaths(int& numPaths, int& steps, double& maturity);
        boost::numeric::ublas::vector<double> getTimeGrid() { return timeGrid_; };
        boost::numeric::ublas::matrix<double> getMonteCarloPaths() { return paths_; };
};

class BlackKarasinskiModel::BlackKarasinskiProcess: public BlackKarasinski::ShortRateDynamics {
    public:
        BlackKarasinskiProcess(const Parameter& fitting, Real alpha, Real sigma) : ShortRateDynamics(boost::shared_ptr<StochasticProcess1D>(
                                 new OrnsteinUhlenbeckProcess(alpha, sigma))),
          fitting_(fitting) {}

        Real variable(Time t, Rate r) const {
            return std::log(r) - fitting_(t);
        }

        Real shortRate(Time t, Real x) const {
            return std::exp(x + fitting_(t));
        }
    private:
        Parameter fitting_;
};

class BlackKarasinskiModel::BlackKarasinskiHelper {
    public:
        BlackKarasinskiHelper(Size i, Real xMin, Real dx, Real discountBondPrice, const boost::shared_ptr<BlackKarasinski::ShortRateTree>& tree) : size_(tree->size(i)),
          dt_(tree->timeGrid().dt(i)),
          xMin_(xMin), dx_(dx),
          statePrices_(tree->statePrices(i)),
          discountBondPrice_(discountBondPrice) {}

        Real operator()(Real theta) const {
            Real value = discountBondPrice_;
            Real x = xMin_;
            for (Size j=0; j<size_; j++) {
                Real discount = std::exp(-std::exp(theta+x)*dt_);
                value -= statePrices_[j]*discount;
                x += dx_;
            }
            return value;
        }

    private:
        Size size_;
        Time dt_;
        Real xMin_, dx_;
        const Array& statePrices_;
        Real discountBondPrice_;
};

#endif