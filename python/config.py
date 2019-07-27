import QuantLib as ql

# configuration file to set market convention on currency and instrument level

class Optimizer:
    def __init__(self, type):
        if type == 'LM':
            epsfcn, xtol, gtol = 1e-8, 1e-8, 1e-8
            self.optimizer = ql.LevenbergMarquardt(epsfcn, xtol, gtol)
        maxIter, maxStationary, rootEps, funcEps, gradNormEps = 10000, 100, 1e-6, 1e-8, 1e-8
        self.endCriterion = ql.EndCriteria(maxIter, maxStationary, rootEps, funcEps, gradNormEps)

    def get_optimizer(self):
        return self.optimizer, self.endCriterion


class CountryLevelInstruments:
    def __init__(self, ccy):
        self.ccy = ccy
        if self.ccy == 'USD':
            self.calendar = ql.UnitedStates()
            self.dayCounter = ql.Actual360()
            self.depo = {'settlementDays': 2,
                         'dayCounterConvention': ql.Actual360(),
                         'businessDayConvention': ql.ModifiedFollowing,
                         'EOM': False
                         }
            self.futures = {'settlementDays': 2,
                            'months': 3,
                            'dayCounterConvention': ql.Actual360(),
                            'businessDayConvention': ql.ModifiedFollowing,
                            'EOM': True
                            }
            self.swap = {'settlementDays': 2,
                         'fixedLegFrequency': ql.Semiannual,
                         'fixedLegTenor': ql.Period(6, ql.Months),
                         'fixedLegAdjustment': ql.Unadjusted,
                         'fixedLegDayCounter': ql.Actual360(),
                         'floatingLegFrequency': ql.Quarterly,
                         'floatingLegTenor': ql.Period(3, ql.Months),
                         'floatingLegAdjustment': ql.ModifiedFollowing,
                         'index': ql.USDLibor(ql.Period(3, ql.Months)),
                         'EOM': False
                         }

            self.swaption = {'fittable': [('1Y', '1Y'), ('1Y', '2Y'), ('1Y', '5Y'), ('1Y', '10Y'), ('1Y', '20Y'),
                                          ('2Y', '1Y'), ('2Y', '2Y'), ('2Y', '5Y'), ('2Y', '10Y'), ('2Y', '20Y'),
                                          ('5Y', '1Y'), ('5Y', '2Y'), ('5Y', '5Y'), ('5Y', '10Y'), ('5Y', '20Y'),
                                          ('10Y', '1Y'), ('10Y', '2Y'), ('10Y', '5Y'), ('10Y', '10Y'), ('10Y', '15Y'),
                                          ('20Y', '1Y'), ('20Y', '2Y'), ('20Y', '5Y'), ('15Y', '10Y')]}