#include <ql/quantlib.hpp>
#define BOOST_PENDING_INTEGER_LOG2_HPP
#include <boost/integer/integer_log2.hpp>

using namespace QuantLib;

boost::shared_ptr<IborIndex> makeIndex(std::vector<Date> dates, std::vector<Rate> rates) {
    DayCounter dayCounter = Actual360();
    RelinkableHandle<YieldTermStructure> termStructure;
    boost::shared_ptr<IborIndex> index(new USDLibor(3*Months, termStructure));

    Date todaysDate = index->fixingCalendar().adjust(Date(28,May,2019));
    Settings::instance().evaluationDate() = todaysDate;
    dates[0] = index->fixingCalendar().advance(todaysDate, index->fixingDays(), Days);
    termStructure.linkTo(boost::shared_ptr<YieldTermStructure>(new ZeroCurve(dates, rates, dayCounter)));

    return index;
}


boost::shared_ptr<IborIndex> makeIndex() {
    Calendar cal = UnitedStates();
    std::vector<Date> dates;
    std::vector<Rate> rates;
    Date todaysDate = Date(30,May,2019);
    dates.push_back(cal.advance(todaysDate, 1*Days)); rates.push_back(0.023);
    dates.push_back(cal.advance(todaysDate, 2*Days)); rates.push_back(0.024);
    dates.push_back(cal.advance(todaysDate, 1*Months)); rates.push_back(0.0236);
    dates.push_back(cal.advance(todaysDate, 3*Months)); rates.push_back(0.0247);
    dates.push_back(cal.advance(todaysDate, 3*Years)); rates.push_back(0.02062);
    dates.push_back(cal.advance(todaysDate, 4*Years)); rates.push_back(0.02044);
    dates.push_back(cal.advance(todaysDate, 5*Years)); rates.push_back(0.02051);
    dates.push_back(cal.advance(todaysDate, 6*Years)); rates.push_back(0.02074);
    dates.push_back(cal.advance(todaysDate, 7*Years)); rates.push_back(0.02102);
    dates.push_back(cal.advance(todaysDate, 8*Years)); rates.push_back(0.02133);
    dates.push_back(cal.advance(todaysDate, 9*Years)); rates.push_back(0.02165);
    dates.push_back(cal.advance(todaysDate, 10*Years)); rates.push_back(0.02198);
    dates.push_back(cal.advance(todaysDate, 15*Years)); rates.push_back(0.02315);
    dates.push_back(cal.advance(todaysDate, 20*Years)); rates.push_back(0.02368);
    dates.push_back(cal.advance(todaysDate, 25*Years)); rates.push_back(0.02384);
    dates.push_back(cal.advance(todaysDate, 30*Years)); rates.push_back(0.0239);
    return makeIndex(dates, rates);
}

boost::shared_ptr<OptionletVolatilityStructure> makeCapVolCurve(const Date& todaysDate) {
    Volatility vols[] = {14.40, 17.15, 16.81, 16.64, 16.17, 15.78, 15.40, 15.21, 14.86};

    std::vector<Date> dates;
    std::vector<Volatility> capletVols;
    boost::shared_ptr<LiborForwardModelProcess> process(new LiborForwardModelProcess(10, makeIndex()));

    for (Size i=0; i < 9; ++i) {
        capletVols.push_back(vols[i]/100);
        dates.push_back(process->fixingDates()[i+1]);
    }

    return boost::shared_ptr<CapletVarianceCurve>(new CapletVarianceCurve(todaysDate, dates, capletVols, Actual360()));
}

void testSimpleCovarianceModels() {
    std::cout << "Testing simple covariance models..." << std::endl;

    SavedSettings backup;

    const Size size = 10;
    const Real tolerance = 1e-14;
    Size i;

    boost::shared_ptr<LmCorrelationModel> corrModel(new LmExponentialCorrelationModel(size, 0.1));

    Matrix recon = corrModel->correlation(0.0) - corrModel->pseudoSqrt(0.0)*transpose(corrModel->pseudoSqrt(0.0));

    for (i=0; i<size; ++i) {
        for (Size j=0; j<size; ++j) {
            if (std::fabs(recon[i][j]) > tolerance)
                std::cout << "Failed to reproduce correlation matrix" << "\n    calculated: " << recon[i][j] << "\n    expected:   " << std::endl;
        }
    }

    std::vector<Time> fixingTimes(size);
    for (i=0; i<size; ++i) {
        fixingTimes[i] = 0.5*i;
    }

    const Real a=0.2;
    const Real b=0.1;
    const Real c=2.1;
    const Real d=0.3;

    boost::shared_ptr<LmVolatilityModel> volaModel(new LmLinearExponentialVolatilityModel(fixingTimes, a, b, c, d));

    boost::shared_ptr<LfmCovarianceProxy> covarProxy(new LfmCovarianceProxy(volaModel, corrModel));

    boost::shared_ptr<LiborForwardModelProcess> process(new LiborForwardModelProcess(size, makeIndex()));

    boost::shared_ptr<LiborForwardModel> liborModel(new LiborForwardModel(process, volaModel, corrModel));

    for (Real t=0; t<4.6; t+=0.31) {
        recon = covarProxy->covariance(t) - covarProxy->diffusion(t)*transpose(covarProxy->diffusion(t));

        for (Size i=0; i<size; ++i) {
            for (Size j=0; j<size; ++j) {
                if (std::fabs(recon[i][j]) > tolerance)
                    std::cout << "Failed to reproduce correlation matrix" << "\n    calculated: " << recon[i][j] << "\n    expected:   " << 0 << std::endl;
            }
        }

        Array volatility = volaModel->volatility(t);

        for (Size k=0; k<size; ++k) {
            Real expected = 0;
            if (k>2*t) {
                const Real T = fixingTimes[k];
                expected=(a*(T-t)+d)*std::exp(-b*(T-t)) + c;
            }

            if (std::fabs(expected - volatility[k]) > tolerance)
                std:: cout << "Failed to reproduce volatities" << "\n    calculated: " << volatility[k] << "\n    expected:   " << expected << std::endl;
        }
    }
}


void testCapletPricing() {
    std::cout << "Testing caplet pricing..." << std::endl;

    SavedSettings backup;

    const Size size = 10;
    #if defined(QL_USE_INDEXED_COUPON)
    const Real tolerance = 1e-5;
    #else
    const Real tolerance = 1e-12;
    #endif

    boost::shared_ptr<IborIndex> index = makeIndex();
    boost::shared_ptr<LiborForwardModelProcess> process(new LiborForwardModelProcess(size, index));

    // set-up pricing engine
    const boost::shared_ptr<OptionletVolatilityStructure> capVolCurve = makeCapVolCurve(Settings::instance().evaluationDate());

    Array variances = LfmHullWhiteParameterization(process, capVolCurve).covariance(0.0).diagonal();

    boost::shared_ptr<LmVolatilityModel> volaModel(new LmFixedVolatilityModel(Sqrt(variances), process->fixingTimes()));

    boost::shared_ptr<LmCorrelationModel> corrModel(new LmExponentialCorrelationModel(size, 0.3));

    boost::shared_ptr<AffineModel> model(new LiborForwardModel(process, volaModel, corrModel));

    const Handle<YieldTermStructure> termStructure = process->index()->forwardingTermStructure();

    boost::shared_ptr<AnalyticCapFloorEngine> engine1(new AnalyticCapFloorEngine(model, termStructure));

    boost::shared_ptr<Cap> cap1(new Cap(process->cashFlows(), std::vector<Rate>(size, 0.04)));
    cap1->setPricingEngine(engine1);

    const Real expected = 0.015853935178;
    const Real calculated = cap1->NPV();

    if (std::fabs(expected - calculated) > tolerance)
        std::cout << "Failed to reproduce npv" << "\n    calculated: " << calculated << "\n    expected:   " << expected;
}

void testCalibration() {
    std::cout << "Testing calibration of a Libor forward model..." << std::endl;

    SavedSettings backup;

    const Size size = 40;
    const Real tolerance = 8e-3;

    boost::shared_ptr<IborIndex> index = makeIndex();
    std::cout << "index made." << std::endl;
    boost::shared_ptr<LiborForwardModelProcess> process(new LiborForwardModelProcess(size, index));
    Handle<YieldTermStructure> termStructure = index->forwardingTermStructure();

    // set-up the model
    boost::shared_ptr<LmVolatilityModel> volaModel(new LmExtLinearExponentialVolModel(process->fixingTimes(), 0.5,0.6,0.1,0.1));

    boost::shared_ptr<LmCorrelationModel> corrModel(new LmLinearExponentialCorrelationModel(size, 0.5, 0.8));

    boost::shared_ptr<LiborForwardModel> model(new LiborForwardModel(process, volaModel, corrModel));

    Size swapVolIndex = 0;
    DayCounter dayCounter=index->forwardingTermStructure()->dayCounter();

    // set-up calibration helper
    std::vector<boost::shared_ptr<CalibrationHelper> > calibrationHelper;

    Period maturity = 6*Months, len = 1*Years;
    Handle<Quote> swaptionVol6x12(boost::shared_ptr<Quote>(new SimpleQuote(0.2838)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper6x12(new SwaptionHelper(maturity, len, swaptionVol6x12, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper6x12->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper6x12);
    maturity = 6*Months, len = 2*Years;
    Handle<Quote> swaptionVol6x24(boost::shared_ptr<Quote>(new SimpleQuote(0.3334)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper6x24(new SwaptionHelper(maturity, len, swaptionVol6x24, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper6x24->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper6x24);
    maturity = 6*Months, len = 3*Years;
    Handle<Quote> swaptionVol6x36(boost::shared_ptr<Quote>(new SimpleQuote(0.3372)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper6x36(new SwaptionHelper(maturity, len, swaptionVol6x36, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper6x36->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper6x36);
    maturity = 6*Months, len = 4*Years;
    Handle<Quote> swaptionVol6x48(boost::shared_ptr<Quote>(new SimpleQuote(0.3271)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper6x48(new SwaptionHelper(maturity, len, swaptionVol6x48, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper6x48->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper6x48);
    maturity = 6*Months, len = 5*Years;
    Handle<Quote> swaptionVol6x60(boost::shared_ptr<Quote>(new SimpleQuote(0.3159)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper6x60(new SwaptionHelper(maturity, len, swaptionVol6x60, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper6x60->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    maturity = 6*Months, len = 6*Years;
    Handle<Quote> swaptionVol6x72(boost::shared_ptr<Quote>(new SimpleQuote(0.3061)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper6x72(new SwaptionHelper(maturity, len, swaptionVol6x72, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper6x72->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper6x72);
    maturity = 6*Months, len = 7*Years;
    Handle<Quote> swaptionVol6x84(boost::shared_ptr<Quote>(new SimpleQuote(0.2965)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper6x84(new SwaptionHelper(maturity, len, swaptionVol6x84, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper6x84->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper6x84);
    maturity = 6*Months, len = 8*Years;
    Handle<Quote> swaptionVol6x96(boost::shared_ptr<Quote>(new SimpleQuote(0.2874)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper6x96(new SwaptionHelper(maturity, len, swaptionVol6x96, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper6x96->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper6x96);
    maturity = 6*Months, len = 9*Years;
    Handle<Quote> swaptionVol6x108(boost::shared_ptr<Quote>(new SimpleQuote(0.2788)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper6x108(new SwaptionHelper(maturity, len, swaptionVol6x108, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper6x108->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper6x108);
    maturity = 1*Years, len = 1*Years;
    Handle<Quote> swaptionVol12x12(boost::shared_ptr<Quote>(new SimpleQuote(0.3519)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper12x12(new SwaptionHelper(maturity, len, swaptionVol12x12, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper12x12->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper12x12);
    maturity = 1*Years, len = 2*Years;
    Handle<Quote> swaptionVol12x24(boost::shared_ptr<Quote>(new SimpleQuote(0.3627)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper12x24(new SwaptionHelper(maturity, len, swaptionVol12x24, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper12x24->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper12x24);
    maturity = 1*Years, len = 3*Years;
    Handle<Quote> swaptionVol12x36(boost::shared_ptr<Quote>(new SimpleQuote(0.3501)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper12x36(new SwaptionHelper(maturity, len, swaptionVol12x36, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper12x36->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper12x36);
    maturity = 1*Years, len = 4*Years;
    Handle<Quote> swaptionVol12x48(boost::shared_ptr<Quote>(new SimpleQuote(0.3337)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper12x48(new SwaptionHelper(maturity, len, swaptionVol12x48, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper12x48->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper12x48);
    maturity = 1*Years, len = 5*Years;
    Handle<Quote> swaptionVol12x60(boost::shared_ptr<Quote>(new SimpleQuote(0.3195)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper12x60(new SwaptionHelper(maturity, len, swaptionVol12x60, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper12x60->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper12x60);
    maturity = 1*Years, len = 6*Years;
    Handle<Quote> swaptionVol12x72(boost::shared_ptr<Quote>(new SimpleQuote(0.3083)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper12x72(new SwaptionHelper(maturity, len, swaptionVol12x72, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper12x72->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper12x72);
    maturity = 1*Years, len = 7*Years;
    Handle<Quote> swaptionVol12x84(boost::shared_ptr<Quote>(new SimpleQuote(0.2978)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper12x84(new SwaptionHelper(maturity, len, swaptionVol12x84, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper12x84->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper12x84);
    maturity = 1*Years, len = 8*Years;
    Handle<Quote> swaptionVol12x96(boost::shared_ptr<Quote>(new SimpleQuote(0.2880)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper12x96(new SwaptionHelper(maturity, len, swaptionVol12x96, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper12x96->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper12x96);
    maturity = 1*Years, len = 9*Years;
    Handle<Quote> swaptionVol12x108(boost::shared_ptr<Quote>(new SimpleQuote(0.2789)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper12x108(new SwaptionHelper(maturity, len, swaptionVol12x108, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper12x108->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper12x108);
    maturity = 2*Years, len = 1*Years;
    Handle<Quote> swaptionVol24x12(boost::shared_ptr<Quote>(new SimpleQuote(0.3723)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper24x12(new SwaptionHelper(maturity, len, swaptionVol24x12, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper24x12->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper24x12);
    maturity = 2*Years, len = 2*Years;
    Handle<Quote> swaptionVol24x24(boost::shared_ptr<Quote>(new SimpleQuote(0.3575)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper24x24(new SwaptionHelper(maturity, len, swaptionVol24x24, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper24x24->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper24x24);
    maturity = 2*Years, len = 3*Years;
    Handle<Quote> swaptionVol24x36(boost::shared_ptr<Quote>(new SimpleQuote(0.3414)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper24x36(new SwaptionHelper(maturity, len, swaptionVol24x36, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper24x36->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper24x36);
    maturity = 2*Years, len = 4*Years;
    Handle<Quote> swaptionVol24x48(boost::shared_ptr<Quote>(new SimpleQuote(0.3277)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper24x48(new SwaptionHelper(maturity, len, swaptionVol24x48, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper24x48->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper24x48);
    maturity = 2*Years, len = 5*Years;
    Handle<Quote> swaptionVol24x60(boost::shared_ptr<Quote>(new SimpleQuote(0.3110)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper24x60(new SwaptionHelper(maturity, len, swaptionVol24x60, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper24x60->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper24x60);
    maturity = 2*Years, len = 6*Years;
    Handle<Quote> swaptionVol24x72(boost::shared_ptr<Quote>(new SimpleQuote(0.3007)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper24x72(new SwaptionHelper(maturity, len, swaptionVol24x72, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper24x72->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper24x72);
    maturity = 2*Years, len = 7*Years;
    Handle<Quote> swaptionVol24x84(boost::shared_ptr<Quote>(new SimpleQuote(0.2913)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper24x84(new SwaptionHelper(maturity, len, swaptionVol24x84, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper24x84->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper24x84);
    maturity = 2*Years, len = 8*Years;
    Handle<Quote> swaptionVol24x96(boost::shared_ptr<Quote>(new SimpleQuote(0.2825)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper24x96(new SwaptionHelper(maturity, len, swaptionVol24x96, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper24x96->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper24x96);
    maturity = 3*Years, len = 1*Years;
    Handle<Quote> swaptionVol36x12(boost::shared_ptr<Quote>(new SimpleQuote(0.3568)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper36x12(new SwaptionHelper(maturity, len, swaptionVol36x12, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper36x12->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper36x12);
    maturity = 3*Years, len = 2*Years;
    Handle<Quote> swaptionVol36x24(boost::shared_ptr<Quote>(new SimpleQuote(0.3396)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper36x24(new SwaptionHelper(maturity, len, swaptionVol36x24, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper36x24->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper36x24);
    maturity = 3*Years, len = 3*Years;
    Handle<Quote> swaptionVol36x36(boost::shared_ptr<Quote>(new SimpleQuote(0.3261)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper36x36(new SwaptionHelper(maturity, len, swaptionVol36x36, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper36x36->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper36x36);
    maturity = 3*Years, len = 4*Years;
    Handle<Quote> swaptionVol36x48(boost::shared_ptr<Quote>(new SimpleQuote(0.3121)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper36x48(new SwaptionHelper(maturity, len, swaptionVol36x48, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper36x48->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper36x48);
    maturity = 3*Years, len = 5*Years;
    Handle<Quote> swaptionVol36x60(boost::shared_ptr<Quote>(new SimpleQuote(0.2992)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper36x60(new SwaptionHelper(maturity, len, swaptionVol36x60, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper36x60->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper36x60);
    maturity = 4*Years, len = 1*Years;
    Handle<Quote> swaptionVol48x12(boost::shared_ptr<Quote>(new SimpleQuote(0.3347)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper48x12(new SwaptionHelper(maturity, len, swaptionVol48x12, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper48x12->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper48x12);
    maturity = 4*Years, len = 2*Years;
    Handle<Quote> swaptionVol48x24(boost::shared_ptr<Quote>(new SimpleQuote(0.3222)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper48x24(new SwaptionHelper(maturity, len, swaptionVol48x24, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper48x24->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper48x24);
    maturity = 4*Years, len = 3*Years;
    Handle<Quote> swaptionVol48x36(boost::shared_ptr<Quote>(new SimpleQuote(0.3107)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper48x36(new SwaptionHelper(maturity, len, swaptionVol48x36, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper48x36->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper48x36);
    maturity = 4*Years, len = 4*Years;
    Handle<Quote> swaptionVol48x48(boost::shared_ptr<Quote>(new SimpleQuote(0.3010)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper48x48(new SwaptionHelper(maturity, len, swaptionVol48x48, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper48x48->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper48x48);
    maturity = 4*Years, len = 5*Years;
    Handle<Quote> swaptionVol48x60(boost::shared_ptr<Quote>(new SimpleQuote(0.2887)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper48x60(new SwaptionHelper(maturity, len, swaptionVol48x60, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper48x60->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper48x60);
    maturity = 5*Years, len = 1*Years;
    Handle<Quote> swaptionVol60x12(boost::shared_ptr<Quote>(new SimpleQuote(0.3210)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper60x12(new SwaptionHelper(maturity, len, swaptionVol60x12, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper60x12->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper60x12);
    maturity = 5*Years, len = 2*Years;
    Handle<Quote> swaptionVol60x24(boost::shared_ptr<Quote>(new SimpleQuote(0.3122)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper60x24(new SwaptionHelper(maturity, len, swaptionVol60x24, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper60x24->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper60x24);
    maturity = 5*Years, len = 3*Years;
    Handle<Quote> swaptionVol60x36(boost::shared_ptr<Quote>(new SimpleQuote(0.3017)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper60x36(new SwaptionHelper(maturity, len, swaptionVol60x36, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper60x36->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper60x36);
    maturity = 5*Years, len = 4*Years;
    Handle<Quote> swaptionVol60x48(boost::shared_ptr<Quote>(new SimpleQuote(0.2913)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper60x48(new SwaptionHelper(maturity, len, swaptionVol60x48, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper60x48->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper60x48);
    maturity = 5*Years, len = 5*Years;
    Handle<Quote> swaptionVol60x60(boost::shared_ptr<Quote>(new SimpleQuote(0.2820)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper60x60(new SwaptionHelper(maturity, len, swaptionVol60x60, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper60x60->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper60x60);
maturity = 6*Years, len = 1*Years;
    Handle<Quote> swaptionVol72x12(boost::shared_ptr<Quote>(new SimpleQuote(0.3046)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper72x12(new SwaptionHelper(maturity, len, swaptionVol72x12, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper72x12->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper72x12);
    maturity = 6*Years, len = 2*Years;
    Handle<Quote> swaptionVol72x24(boost::shared_ptr<Quote>(new SimpleQuote(0.2958)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper72x24(new SwaptionHelper(maturity, len, swaptionVol72x24, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper72x24->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper72x24);
    maturity = 6*Years, len = 3*Years;
    Handle<Quote> swaptionVol72x36(boost::shared_ptr<Quote>(new SimpleQuote(0.2876)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper72x36(new SwaptionHelper(maturity, len, swaptionVol72x36, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper72x36->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper72x36);
    maturity = 6*Years, len = 4*Years;
    Handle<Quote> swaptionVol72x48(boost::shared_ptr<Quote>(new SimpleQuote(0.2797)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper72x48(new SwaptionHelper(maturity, len, swaptionVol72x48, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper72x48->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper72x48);
    maturity = 7*Years, len = 1*Years;
    Handle<Quote> swaptionVol84x12(boost::shared_ptr<Quote>(new SimpleQuote(0.2902)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper84x12(new SwaptionHelper(maturity, len, swaptionVol84x12, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper84x12->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper84x12);
    maturity = 7*Years, len = 2*Years;
    Handle<Quote> swaptionVol84x24(boost::shared_ptr<Quote>(new SimpleQuote(0.2813)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper84x24(new SwaptionHelper(maturity, len, swaptionVol84x24, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper84x24->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper84x24);
    maturity = 7*Years, len = 3*Years;
    Handle<Quote> swaptionVol84x36(boost::shared_ptr<Quote>(new SimpleQuote(0.2750)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper84x36(new SwaptionHelper(maturity, len, swaptionVol84x36, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper84x36->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper84x36);
    maturity = 8*Years, len = 1*Years;
    Handle<Quote> swaptionVol96x12(boost::shared_ptr<Quote>(new SimpleQuote(0.2776)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper96x12(new SwaptionHelper(maturity, len, swaptionVol96x12, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper96x12->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper96x12);
    maturity = 8*Years, len = 2*Years;
    Handle<Quote> swaptionVol96x24(boost::shared_ptr<Quote>(new SimpleQuote(0.2707)));
    boost::shared_ptr<CalibrationHelper> swaptionHelper96x24(new SwaptionHelper(maturity, len, swaptionVol96x24, index, index->tenor(), dayCounter, index->dayCounter(), termStructure, CalibrationHelper::ImpliedVolError));
    swaptionHelper96x24->setPricingEngine(boost::shared_ptr<PricingEngine>(new LfmSwaptionEngine(model,termStructure)));
    calibrationHelper.push_back(swaptionHelper96x24);


    LevenbergMarquardt om(1e-6, 1e-6, 1e-6);
    model->calibrate(calibrationHelper, om, EndCriteria(2000, 100, 1e-6, 1e-6, 1e-6));
    std::cout << model.get()->params() << std::endl;
    // measure the calibration error
    Real calculated = 0.0;
    for (Size i=0; i<calibrationHelper.size(); ++i) {
        Real diff = calibrationHelper[i]->calibrationError();
        std::cout << "Calibration Instrument " << i << ", calibration error: " << diff << std::endl;
        calculated += diff*diff;
    }

    //if (std::sqrt(calculated) > tolerance)
    //    std::cout << "Failed to calibrate libor forward model" << "\n    calculated diff: " << std::sqrt(calculated) << "\n    expected : smaller than  " << tolerance << std::endl;
}

void BiInterp(){
    Date todaysDate = UnitedStates().adjust(Date(28,May,2019));
    std::vector<Period> expiry, tenor;
    std::vector<Real> expiries, tenors;
    expiry.push_back(1*Months); expiry.push_back(3*Months); expiry.push_back(6*Months); expiry.push_back(1*Years); expiry.push_back(2*Years); expiry.push_back(3*Years); expiry.push_back(4*Years); expiry.push_back(5*Years); expiry.push_back(6*Years);  expiry.push_back(7*Years); expiry.push_back(8*Years); expiry.push_back(9*Years); expiry.push_back(10*Years); expiry.push_back(15*Years); expiry.push_back(20*Years);  expiry.push_back(25*Years); expiry.push_back(30*Years);
    tenor.push_back(1*Years); tenor.push_back(2*Years); tenor.push_back(3*Years); tenor.push_back(4*Years); tenor.push_back(5*Years); tenor.push_back(6*Years);  tenor.push_back(7*Years); tenor.push_back(8*Years); tenor.push_back(9*Years); tenor.push_back(10*Years); tenor.push_back(15*Years); tenor.push_back(20*Years);  tenor.push_back(25*Years); tenor.push_back(30*Years); 
    Volatility swaptionVols[] = {0.1985, 0.2895, 0.3031, 0.3093, 0.3097, 0.3026, 0.2950, 0.2834, 0.2721, 0.2610, 0.2390, 0.2207, 0.2125, 0.2050,
                                 0.2224, 0.3037, 0.3086, 0.3108, 0.3090, 0.2991, 0.2892, 0.2795, 0.2704, 0.2615, 0.2365, 0.2223, 0.2169, 0.2097,
                                 0.2838, 0.3334, 0.3372, 0.3271, 0.3159, 0.3061, 0.2965, 0.2874, 0.2788, 0.2705, 0.2439, 0.2287, 0.2228, 0.2167,
                                 0.3519, 0.3627, 0.3501, 0.3337, 0.3195, 0.3083, 0.2978, 0.2880, 0.2789, 0.2704, 0.2441, 0.2295, 0.2233, 0.2190,
                                 0.3723, 0.3575, 0.3414, 0.3277, 0.3110, 0.3007, 0.2913, 0.2825, 0.2742, 0.2669, 0.2419, 0.2284, 0.2243, 0.2203,
                                 0.3568, 0.3396, 0.3261, 0.3121, 0.2992, 0.2903, 0.2822, 0.2742, 0.2673, 0.2610, 0.2387, 0.2258, 0.2229, 0.2198,
                                 0.3347, 0.3222, 0.3107, 0.3010, 0.2887, 0.2810, 0.2737, 0.2671, 0.2610, 0.2557, 0.2354, 0.2241, 0.2215, 0.2192,
                                 0.3210, 0.3122, 0.3017, 0.2913, 0.2820, 0.2747, 0.2686, 0.2627, 0.2576, 0.2529, 0.2341, 0.2237, 0.2217, 0.2201,
                                 0.3046, 0.2958, 0.2876, 0.2797, 0.2719, 0.2662, 0.2610, 0.2562, 0.2518, 0.2478, 0.2310, 0.2217, 0.2199, 0.2191,
                                 0.2902, 0.2813, 0.2750, 0.2688, 0.2640, 0.2591, 0.2550, 0.2510, 0.2473, 0.2437, 0.2285, 0.2201, 0.2187, 0.2185,
                                 0.2776, 0.2707, 0.2641, 0.2603, 0.2568, 0.2531, 0.2496, 0.2461, 0.2428, 0.2397, 0.2257, 0.2184, 0.2171, 0.2178, 
                                 0.2676, 0.2606, 0.2565, 0.2536, 0.2514, 0.2481, 0.2452, 0.2422, 0.2394, 0.2367, 0.2238, 0.2171, 0.2163, 0.2175, 
                                 0.2561, 0.2531, 0.2500, 0.2485, 0.2471, 0.2444, 0.2418, 0.2393, 0.2368, 0.2344, 0.2223, 0.2162, 0.2154, 0.2176,
                                 0.2434, 0.2403, 0.2392, 0.2385, 0.2380, 0.2365, 0.2351, 0.2335, 0.2318, 0.2300, 0.2208, 0.2170, 0.2193, 0.2226,
                                 0.2378, 0.2365, 0.2353, 0.2345, 0.2333, 0.2318, 0.2302, 0.2283, 0.2262, 0.2243, 0.2183, 0.2194, 0.2231, 0.2250,
                                 0.2362, 0.2357, 0.2339, 0.2324, 0.2310, 0.2295, 0.2281, 0.2267, 0.2256, 0.2247, 0.2262, 0.2299, 0.2316, 0.2337,
                                 0.2274, 0.2280, 0.2278, 0.2273, 0.2271, 0.2280, 0.2291, 0.2302, 0.2313, 0.2322, 0.2374, 0.2383, 0.2402, 0.2431};
    //Interpolation2D interp = FlatExtrapolator2D(boost::make_shared<BilinearInterpolation>(expiries.begin(), expiries.end(), tenors.begin(), tenors.end(), swaptionVols));
    //std::cout << "6MX3M" << interp.
    return;
}