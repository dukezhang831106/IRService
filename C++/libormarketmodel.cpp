#include "libormarketmodel.hpp"

void LiborMarketModel::parse_calendar(std::string& holmasks) {
    std::stringstream holidays(holmasks);
    int cnt = 0;
    while(holidays.good()){
        std::string hol;
        getline(holidays, hol, ',');
        Calendar cal;
        if (hol == "NY") 
            cal = UnitedStates();
        else if (hol == "LON")
            cal = UnitedKingdom();
        else if (hol == "EUR")
            cal = TARGET();
        else if (hol == "TOK")
            cal = Japan();
        else if (hol == "BEI")
            cal = China();
        else if (hol == "SYD")
            cal = Australia();
        else if (hol == "ZUR")
            cal = Switzerland();
        else 
            cal = UnitedStates();
        if (cnt == 0){
            calendar_ = cal;
        }
        else{
            calendar_ = JointCalendar(calendar_, cal);
        }
        cnt++;
    }
}

LiborMarketModel::LiborMarketModel(std::string modelDate, std::string currency = "USD", std::string exchange = "NY"){
    
}



void LiborMarketModel::buildCurve(){
    std::ifstream file("data/configuration.json");
    json output;
    file >> output;
    for(auto & exchange : output[currency_]){
        if (exchange["Exchange"] == exchange_){
            //Calendar cal = parse_calendar(std::string(exchange["Holiday"]));
            std::cout << exchange["Instruments"] << std::endl;
        }
    }
    return;
}

//void LiborMarketModel::calibration(){
//}

void LiborMarketModel::test_connection(){
    pqxx::connection conn("dbname=EikonInstrument user=postgres host=127.0.0.1 port=1111 password=123456");
    pqxx::work txn(conn);

    std::stringstream command;
    command << "SELECT * " <<
        "FROM \"IRInstruments\" " <<
        "WHERE \"QuoteDate\"=Date('" << modelModelDate_ << "') AND \"Currency\"='USD' AND \"InstrumentType\"='Swaption';";

    pqxx::result res = txn.exec(command.str());
    
    for(int i = 0; i < res.size(); i++){
        for(int j = 0; j < res[i].size(); j++){
            std::cout << res[i][j] << ", ";
        }
        std::cout << std::endl;
    }
    return;
}