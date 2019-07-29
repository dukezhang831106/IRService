#include "libormarketmodel.hpp"
//#include <iostream>

int main(){
    std::string currency = "USD", holmask = "NY", modelDate = "2019-05-28";
    LiborMarketModel model(modelDate, currency, holmask);
    //std::string freq = "7D";
    //Period f = model.parse_period(freq);
    //std::cout << f << std::endl;
    //model.buildCurve();
    //model.test_connection();
    return EXIT_SUCCESS;
}