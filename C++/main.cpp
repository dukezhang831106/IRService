#include "libormarketmodel.hpp"
#include "shortratemodel.hpp"

int main(int argc, char** argv){
    if (argc > 2){
        std::string currency = "USD", exchange = "NY", interpolationType = "LogLinear";
        std::string modelType = argv[1], modelDate = argv[2];
        std::cout << "Building " << modelType << " on " << modelDate << std::endl;;
        if (modelType == "HullWhite" || modelType == "HW"){
            HullWhiteModel hw(modelDate, currency, exchange, interpolationType);
            std::vector<CalibrationReport> report = hw.get_calibrationReport();
            for(auto& v : report){
                std::cout << v.instrumentType << " " << v.expiry << "X" << v.tenor << ", expected: " << v.expected << ", implied: " << v.calculated << ", diff: " << v.diff << std::endl;
            }
            int numPaths = 20, timeSteps = 10;
            double maturity = 30.0;
            hw.generateMonteCarloPaths(numPaths, timeSteps, maturity);
            boost::numeric::ublas::vector<double> timeGrids = hw.getTimeGrid();
            boost::numeric::ublas::matrix<double> ratePaths = hw.getMonteCarloPaths();
            std::cout << "Time Grid: " << timeGrids << std::endl;
            std::cout << "Simulated Paths: " << ratePaths << std::endl;
            return EXIT_SUCCESS;
        }
        else if (modelType == "CoxIngersollRoss" || modelType == "CIR"){
            CoxIngersollRossModel cir(modelDate, currency, exchange, interpolationType);
            std::vector<CalibrationReport> report = cir.get_calibrationReport();
            for(auto& v : report){
                std::cout << v.instrumentType << " " << v.expiry << "X" << v.tenor << ", expected: " << v.expected << ", implied: " << v.calculated << ", diff: " << v.diff << std::endl;
            }
            std::vector<Period> maturities = {1*Years, 2*Years, 3*Years, 5*Years, 7*Years, 10*Years, 12*Years, 15*Years, 20*Years};
            ProbabilityTree ptree = cir.getLattice(maturities);
            for(int level = 0; level < ptree.R.size(); level++){
                std::cout << "Level " << level << std::endl;
                std::cout << "forward: " << ptree.Fwd[level] << std::endl;
                if (level < ptree.R.size() - 1){
                    std::cout << "connection: " << ptree.Connect[level] << std::endl;
                    std::cout << "prob: " << ptree.Prob[level] << std::endl;
                }
            }
            return EXIT_SUCCESS;
        }
        else if (modelType == "BlackKarasinski" || modelType == "BK"){
            BlackKarasinskiModel bk(modelDate, currency, exchange, interpolationType);
            std::vector<CalibrationReport> report = bk.getCalibrationReport();
            for(auto& v : report){
                std::cout << v.instrumentType << " " << v.expiry << "X" << v.tenor << ", expected: " << v.expected << ", implied: " << v.calculated << ", diff: " << v.diff << std::endl;
            }
            std::vector<Period> maturities = {1*Years, 2*Years, 3*Years, 5*Years, 7*Years, 10*Years, 12*Years, 15*Years, 20*Years};
            ProbabilityTree ptree = bk.getLattice(maturities);
            for(int level = 0; level < maturities.size(); level++){
                std::cout << "Level: " << level << ", Date: " << ptree.matDates[level] << ", fwd rate: " << ptree.R[level] << std::endl;
                if (level > 0){
                    std::cout << std::endl;
                    std::cout << "prob " << level-1 << " to " << level << ": " << ptree.Prob[level-1] << std::endl;
                }
            }
            return EXIT_SUCCESS;
        }
        else if (modelType == "LMM" || modelType == "BGM"){
            int numPaths = 3, timeSteps = 10;
            double maturity = 30.0;
            LiborMarketModel lmm(modelDate, currency, exchange, interpolationType);
            std::vector<CalibrationReport> report = lmm.get_calibrationReport();
            for(auto& v : report){
                std::cout << v.instrumentType << " " << v.expiry << "X" << v.tenor << ", expected: " << v.expected << ", implied: " << v.calculated << ", diff: " << v.diff << std::endl;
            }
            lmm.generateMonteCarloPaths(numPaths, timeSteps, maturity);
            boost::numeric::ublas::vector<double> timeGrids = lmm.getTimeGrid();
            boost::numeric::ublas::matrix<double> ratePaths = lmm.getMonteCarloPaths();
            std::cout << "Time Grid: " << timeGrids << std::endl;
            std::cout << "Simulated Paths: " << ratePaths << std::endl;
            return EXIT_SUCCESS;
        }
    }
    return EXIT_FAILURE;
}