#include "libormarketmodel.hpp"
//#include <iostream>

int main(){
    std::string currency = "USD", holmask = "NY", modelDate = "2019-05-28";
    LiborMarketModel model(modelDate, currency, holmask);
    model.buildModel();
    //model.test_connection();
    return EXIT_SUCCESS;
}