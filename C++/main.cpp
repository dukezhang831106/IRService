#include "libormarketmodel.hpp"
#include "shortratemodel.hpp"

int main(){
    std::string currency = "USD", exchange = "NY", modelDate = "2019-05-28", interpolationType = "LogLinear";
    
    HullWhiteModel hw(modelDate, currency, exchange, interpolationType);
    //ExtendedCIRModel cir(modelDate, currency, exchange, interpolationType);
    //BlackKarasinskiModel bk(modelDate, currency, exchange, interpolationType);
    //std::string dayCounter = "30/360", interpolationType = "LogLinear";
    //model.buildCurve(dayCounter, interpolationType);
    //std::string freq = "7D";
    //Period f = model.parse_period(freq);
    //std::cout << f << std::endl;
    //model.buildCurve();
    //model.parse_instruments();
    return EXIT_SUCCESS;
}
/*
    boost::shared_ptr<IborIndex> makeIndex(std::vector<Date> dates,
                                           std::vector<Rate> rates) {
        DayCounter dayCounter = Actual360();

        RelinkableHandle<YieldTermStructure> termStructure;

        boost::shared_ptr<IborIndex> index(new Euribor6M(termStructure));

        Date todaysDate =
            index->fixingCalendar().adjust(Date(4,September,2005));
        Settings::instance().evaluationDate() = todaysDate;

        dates[0] = index->fixingCalendar().advance(todaysDate,
                                                   index->fixingDays(), Days);

        termStructure.linkTo(boost::shared_ptr<YieldTermStructure>(
                                    new ZeroCurve(dates, rates, dayCounter)));

        return index;
    }


    boost::shared_ptr<IborIndex> makeIndex() {
        std::vector<Date> dates;
        std::vector<Rate> rates;
        dates.push_back(Date(4,September,2005));
        dates.push_back(Date(4,September,2018));
        rates.push_back(0.039);
        rates.push_back(0.041);

        return makeIndex(dates, rates);
    }



int main(){
    SavedSettings backup;

    const Size size = 14;
    const Real tolerance = 8e-3;

    Volatility capVols[] = {0.145708,0.158465,0.166248,0.168672,
                            0.169007,0.167956,0.166261,0.164239,
                            0.162082,0.159923,0.157781,0.155745,
                            0.153776,0.151950,0.150189,0.148582,
                            0.147034,0.145598,0.144248};

    Volatility swaptionVols[] = {0.170595, 0.166844, 0.158306, 0.147444, 0.136930, 0.126833, 0.118135, 
                                 0.175963, 0.166359, 0.155203, 0.143712, 0.132769, 0.122947, 0.114310, 
                                 0.174455, 0.162265, 0.150539, 0.138734, 0.128215, 0.118470, 0.110540, 
                                 0.169780, 0.156860, 0.144821, 0.133537, 0.123167, 0.114363, 0.106500,
                                 0.164521, 0.151223, 0.139670, 0.128632, 0.119123, 0.110330, 0.103114, 
                                 0.158956, 0.146036, 0.134555, 0.124393, 0.115038, 0.106996, 0.100064};

    boost::shared_ptr<IborIndex> index = makeIndex();
    boost::shared_ptr<LiborForwardModelProcess> process(
        new LiborForwardModelProcess(size, index));
    Handle<YieldTermStructure> termStructure = index->forwardingTermStructure();

    // set-up the model
    boost::shared_ptr<LmVolatilityModel> volaModel(
                    new LmExtLinearExponentialVolModel(process->fixingTimes(),
                                                       0.5,0.6,0.1,0.1));

    boost::shared_ptr<LmCorrelationModel> corrModel(
                     new LmLinearExponentialCorrelationModel(size, 0.5, 0.8));

    boost::shared_ptr<LiborForwardModel> model(
                        new LiborForwardModel(process, volaModel, corrModel));

    Size swapVolIndex = 0;
    DayCounter dayCounter=index->forwardingTermStructure()->dayCounter();

    // set-up calibration helper
    std::vector<boost::shared_ptr<CalibrationHelper> > calibrationHelper;

    Size i;
    for (i=2; i < size; ++i) {
        const Period maturity = i*index->tenor();
        Handle<Quote> capVol(
            boost::shared_ptr<Quote>(new SimpleQuote(capVols[i-2])));
        std::cout << "CapFloor(" << maturity << "): " << capVols[i-2] << std::endl;
        boost::shared_ptr<CalibrationHelper> caphelper(
            new CapHelper(maturity, capVol, index, Annual,
                          index->dayCounter(), true, termStructure,
                          CalibrationHelper::ImpliedVolError));

        caphelper->setPricingEngine(boost::shared_ptr<PricingEngine>(
                           new AnalyticCapFloorEngine(model, termStructure)));

        calibrationHelper.push_back(caphelper);

        if (i<= size/2) {
            // add a few swaptions to test swaption calibration as well
            for (Size j=1; j <= size/2; ++j) {
                const Period len = j*index->tenor();
                int q = swapVolIndex++;
                Handle<Quote> swaptionVol(
                    boost::shared_ptr<Quote>(
                        new SimpleQuote(swaptionVols[q])));
                
                std::cout << "Swaptiong(" << maturity << ", " << len << "): " << swaptionVols[q] << std::endl;

                boost::shared_ptr<CalibrationHelper> swaptionHelper(
                    new SwaptionHelper(maturity, len, swaptionVol, index,
                                       index->tenor(), dayCounter,
                                       index->dayCounter(),
                                       termStructure,
                                       CalibrationHelper::ImpliedVolError));

                swaptionHelper->setPricingEngine(
                     boost::shared_ptr<PricingEngine>(
                                 new LfmSwaptionEngine(model,termStructure)));

                calibrationHelper.push_back(swaptionHelper);
            }
        }
    }

    LevenbergMarquardt om(1e-6, 1e-6, 1e-6);
    model->calibrate(calibrationHelper, om, EndCriteria(2000, 100, 1e-6, 1e-6, 1e-6));

    // measure the calibration error
    Real calculated = 0.0;
    for (i=0; i<calibrationHelper.size(); ++i) {
        Real diff = calibrationHelper[i]->calibrationError();
        calculated += diff*diff;
    }

    return EXIT_SUCCESS;
}*/