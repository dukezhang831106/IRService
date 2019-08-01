import QuantLib as ql
from collections import namedtuple
import numpy as np
import math
import psycopg2
import matplotlib.pyplot as plt
from config import *

def generate_bucket_mapping(expiry):
    periods = {'D': ql.Days, 'W': ql.Weeks, 'M': ql.Months, 'Y': ql.Years}
    if ql.IMM_isIMMcode(expiry) == True:
        return ql.IMM_date(expiry)
    elif expiry == 'ON':
        return ql.Period(0, ql.Days)
    elif expiry == 'TN':
        return ql.Period(1, ql.Days)
    elif expiry == 'SN':
        return ql.Period(2, ql.Days)
    return ql.Period(int(expiry[:-1]), periods[expiry[-1]])



class IRCurve:
    def __init__(self, quote_date, ccy):
        self.ccy = ccy
        self.quote_date = quote_date
        self.configuration = CountryLevelInstruments(ccy)
        self.settlement_date = self.configuration.calendar.advance(ql.DateParser_parse(self.quote_date, '%Y-%m-%d'), ql.Period(self.configuration.swap['settlementDays'], ql.Days))
        ql.Settings.instance().evaluationDate = self.settlement_date

    def calibration(self):
        connection = psycopg2.connect(dbname='EikonInstrument', user='postgres', host='127.0.0.1', port='5432',
                                      password='123456')
        cursor = connection.cursor()
        cursor.execute("SELECT * FROM \"IRInstruments\" WHERE \"QuoteDate\"=Date('{}') AND \"Currency\"='{}';".format(self.quote_date, self.ccy))
        instruments = cursor.fetchall()
        connection.close()
        depo, repo, futures, swap, bond = dict(), dict(), dict(), dict(), dict()

        for instrument in instruments:
            if instrument[2] == 'Depo':
                depo[generate_bucket_mapping(instrument[4])] = instrument[6]/100.0
            elif instrument[2] == 'Repo':
                repo[generate_bucket_mapping(instrument[4])] = instrument[6]/100.0
            elif instrument[2] == 'Fut':
                futures[generate_bucket_mapping(instrument[4])] = instrument[6]
            elif instrument[2] == 'Swap':
                swap[generate_bucket_mapping(instrument[4])] = instrument[6]/100.0

        if self.ccy == 'USD':
            depositHelpers = [ql.DepositRateHelper(ql.QuoteHandle(ql.SimpleQuote(depo[key])), key,
                                                   self.configuration.depo['settlementDays'],
                                                   self.configuration.calendar,
                                                   self.configuration.depo['businessDayConvention'],
                                                   self.configuration.depo['EOM'],
                                                   self.configuration.depo['dayCounterConvention']) for key
                              in depo.keys()]

            futuresHelpers = [ql.FuturesRateHelper(ql.QuoteHandle(ql.SimpleQuote(futures[key])), key,
                              self.configuration.futures['months'], self.configuration.calendar,
                              self.configuration.futures['businessDayConvention'],
                              self.configuration.futures['EOM'],
                              self.configuration.futures['dayCounterConvention'],
                              ql.QuoteHandle(ql.SimpleQuote(0.0))) for key in futures.keys()]

            swapHelpers = [ql.SwapRateHelper(swap[key], key, self.configuration.calendar,
                           self.configuration.swap['fixedLegFrequency'],
                       self.configuration.swap['fixedLegAdjustment'],
                       self.configuration.swap['fixedLegDayCounter'],
                       self.configuration.swap['index']) for key in swap.keys()]

        helpers = depositHelpers + futuresHelpers + swapHelpers
        self.discountTermStructure = ql.RelinkableYieldTermStructureHandle()
        self.forecastTermStructure = ql.RelinkableYieldTermStructureHandle()
        for hh in helpers:
            print("{}: {}".format(hh.maturityDate(), hh.quote().value()))
        self.curve = ql.PiecewiseLinearForward(self.settlement_date, helpers, self.configuration.dayCounter)
        self.discountTermStructure.linkTo(self.curve)

    def get_spot_rates(self):
        spots, tenors = list(), list()
        for d in self.curve.dates():
            yrs = self.configuration.dayCounter.yearFraction(self.curve.dates()[0], d)
            compounding = ql.Compounded
            freq = ql.Semiannual
            zero_rate = self.curve.zeroRate(yrs, compounding, freq)
            tenors.append(yrs)
            eq_rate = zero_rate.equivalentRate(self.configuration.dayCounter, compounding, freq, self.curve.dates()[0], d).rate()
            spots.append(100 * eq_rate)
        return tenors, spots

    def get_curve(self):
        return self.curve, self.discountTermStructure

class IRModel:
    def __init__(self, modelType, quote_date, ccy):
        self.modelType = modelType
        self.ccy = ccy
        self.quote_date = quote_date
        self.model_date = ql.DateParser_parse(quote_date, '%Y-%m-%d')
        self.configuration = CountryLevelInstruments(self.ccy)
        self.settlement_date = ql.UnitedStates().advance(self.model_date,
                                                    ql.Period(self.configuration.swap['settlementDays'], ql.Days))
        self.model_name = '{}_{}_{}'.format(self.modelType, self.ccy, self.quote_date)
        ql.Settings.instance().evaluationDate = self.settlement_date
        self.optim = Optimizer('LM') # Levenberg-Marquardt
        self.optimizer = self.optim.optimizer
        self.endCriterion = self.optim.endCriterion

    def update_curve(self, curve, discountTermStructure):
        self.curve = curve
        self.discountTermStructure = discountTermStructure

    def calibration(self):
        index = ql.USDLibor(ql.Period(3, ql.Months), self.discountTermStructure)
        connection = psycopg2.connect(dbname='EikonInstrument', user='postgres', host='127.0.0.1', port='5432',
                                      password='123456')
        cursor = connection.cursor()
        cursor.execute("SELECT * FROM \"IRInstruments\" WHERE \"QuoteDate\"=Date('{}') AND \"Currency\"='{}' AND \"InstrumentType\"='Swaption';".format(
            self.quote_date, self.ccy))
        instruments = cursor.fetchall()
        connection.close()

        CalibrationData = namedtuple("CalibrationData",
                                     "start, length, volatility")

        data = list()
        for instrument in instruments:
            bucket = (instrument[4], instrument[5])
            if bucket in self.configuration.swaption['fittable']:
                data.append(CalibrationData(generate_bucket_mapping(instrument[4]), generate_bucket_mapping(instrument[5]), instrument[6]/100.0))

        if self.modelType == 'HullWhite':
            self.model = ql.HullWhite(ql.YieldTermStructureHandle(self.curve))
            self.engine = ql.JamshidianSwaptionEngine(self.model)
        elif self.modelType == 'BlackKarasinski':
            self.model = ql.BlackKarasinski(ql.YieldTermStructureHandle(self.curve))
            self.engine = ql.TreeSwaptionEngine(self.model, 100)

        elif self.modelType == 'G2++':
            self.model = ql.G2(ql.YieldTermStructureHandle(self.curve))
            self.engine = ql.TreeSwaptionEngine(self.model, 25)

        swaptions = list()
        for d in data:
            vol_handle = ql.QuoteHandle(ql.SimpleQuote(d.volatility))
            helper = ql.SwaptionHelper(d.start, d.length, vol_handle, index, self.configuration.swap['fixedLegTenor'],
                                       self.configuration.swap['fixedLegDayCounter'],
                                       self.configuration.dayCounter, ql.YieldTermStructureHandle(self.curve))
            helper.setPricingEngine(self.engine)
            swaptions.append(helper)

        self.model.calibrate(swaptions, self.optimizer, self.endCriterion)

    def generate_monte_carlo_paths(self, timestep, length, numPaths, lowDiscrepancy=False, brownianBridge=True):
        arr = np.zeros((numPaths, timestep+1))

        if self.modelType == 'HullWhite':
            a, sigma = self.model.params()
            self.process = ql.HullWhiteProcess(ql.YieldTermStructureHandle(self.curve), a, sigma)
        elif self.modelType == 'BlackKarasinski':
            a, sigma = self.model.params()
            #self.process = ql.ExtendedOrnsteinUhlenbeckProcess(a, sigma, )
        elif self.modelType == 'LMM':
            return

        if lowDiscrepancy:
            usg = ql.UniformLowDiscrepancySequenceGenerator(timestep)
            rng = ql.GaussianLowDiscrepancySequenceGenerator(usg)
            seq = ql.GaussianSobolPathGenerator(self.process, length, timestep, rng, brownianBridge)
        else:
            usg = ql.UniformRandomSequenceGenerator(timestep, ql.UniformRandomGenerator())
            rng = ql.GaussianRandomSequenceGenerator(usg)
            seq = ql.GaussianPathGenerator(self.process, length, timestep, rng, brownianBridge)
        maximum, minimum = list(), list()
        for i in range(numPaths):
            sample = seq.next()
            path = sample.value()
            time = [path.time(j) for j in range(len(path))]
            value = [path[j] for j in range(len(path))]
            arr[i, :] = 100.0*np.array(value)
            maximum.append(np.max(arr[i, :]))
            minimum.append(np.min(arr[i, :]))
        return time, arr, np.max(maximum), np.min(minimum)

if __name__ == '__main__':
    quote_date, ccy, model_type = '2019-05-28', 'USD', 'HullWhite'
    curveHandle = IRCurve(quote_date, ccy)
    curveHandle.calibration()
    tenor, spots = curveHandle.get_spot_rates()
    print(tenor)
    print(spots)
    '''curve, termStructure = curveHandle.get_curve()
    tenor, spots = curveHandle.get_spot_rates()
    modelHandle = IRModel(model_type, quote_date, ccy)
    modelHandle.update_curve(curve, termStructure)
    modelHandle.calibration()
    numPaths = 50
    time, paths, maxima, minima = modelHandle.generate_monte_carlo_paths(10, 30, numPaths)
    time = np.array(time)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (12, 4.5))
    ax1.plot(tenor[1:], spots[1:])
    ax1.set_xlabel('Tenor')
    ax1.set_ylabel('Spot Rate')
    ax1.set_title('IR Curve')

    for i in range(numPaths):
        ax2.plot(time, paths[i, :])
    ax2.set_xlabel('Tenor')
    ax2.set_ylabel('Rate')
    ax2.set_title('IR Model Paths')
    plt.show()'''