#ifndef libor_market_model_hpp
#define libor_market_model_hpp

 // avoid pragma warning
#include "basemodel.hpp"

class LiborMarketModel : public BaseModel{
    private:
        std::set<std::string> fittable_;
        std::vector<boost::shared_ptr<CalibrationHelper>> modelCalibrator_;
        boost::shared_ptr<IborIndex>  index_;
        boost::shared_ptr<LiborForwardModel> LMMmodel_;
        // For calibration;
        std::vector<CalibrationReport> calibrationReports_;
        // For Monte Carlo Simulation
        boost::numeric::ublas::vector<double> timeGrid_;
        boost::numeric::ublas::matrix<double> paths_;
        void parse_fitting_bucket(std::string& key, json& instrument);
        void parse_vol_instruments(std::string& key, json& instrument);

    public:
        LiborMarketModel(std::string &modelDate, std::string &currency, std::string &exchange, std::string& interpolationType, bool& hello) : BaseModel(modelDate, currency, exchange, interpolationType) {
            std::string indices = "USD3M";
            index_ = makeIndex(indices);
            const Size size = 120;
            boost::shared_ptr<LiborForwardModelProcess> process(new LiborForwardModelProcess(size, index_));
            // set-up the model
            boost::shared_ptr<LmVolatilityModel> volModel(new LmExtLinearExponentialVolModel(process->fixingTimes(), 0.5,0.6,0.1,0.1));
            boost::shared_ptr<LmCorrelationModel> corrModel(new LmLinearExponentialCorrelationModel(size, 0.5, 0.8));
            LMMmodel_ = boost::shared_ptr<LiborForwardModel>(new LiborForwardModel(process, volModel, corrModel));
        };
        LiborMarketModel(std::string &modelDate, std::string &currency, std::string &exchange, std::string& interpolationType);
        void buildModel();
        std::vector<CalibrationReport> get_calibrationReport() { return calibrationReports_; };
        void generateMonteCarloPaths(int& nodes, int& steps, double& maturity);
        boost::numeric::ublas::vector<double> getTimeGrid() { return timeGrid_; };
        boost::numeric::ublas::matrix<double> getMonteCarloPaths() { return paths_; };
};

#endif