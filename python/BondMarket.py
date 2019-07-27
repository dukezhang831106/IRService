import os
import time
import datetime as dt
import pandas as pd
import requests
import urllib
import urllib3
import json
import QuantLib as ql
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
from concurrent.futures import ThreadPoolExecutor

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

units = {'Week': ql.Weeks, 'Year': ql.Years}
freq = {'None': ql.Period(ql.Weekly), 'Quarterly': ql.Period(ql.Quarterly), 'Semi-Annual': ql.Period(ql.Semiannual),
        'Annual': ql.Period(ql.Annual)}


# bond market
def grab_bond_data(date, market='US_Treasury'):
    data = {
        'priceDate.month': str(date.month),
        'priceDate.day': str(date.day),
        'priceDate.year': str(date.year),
        'submit': 'Show Prices'
    }

    responses = requests.post('https://www.treasurydirect.gov/GA-FI/FedInvest/selectSecurityPriceDate.htm', data=data,
                              verify=False).text.split('<table class="data1">')[1].split('</table>')[0].split(
        '</tr>\r\n\t\r\n\t<tr>')[1:]
    cusips, types, coupons, maturities, calldates, bids, asks, eods = list(), list(), list(), list(), list(), list(), \
                                                                      list(), list()

    r0 = responses[0].split('</td>\r\n\t\t<td>')
    for r in responses:
        response = r.split('</td>\r\n\t\t<td>')
        cusips.append(response[0].split('<td>')[1].split('</td>')[0])
        types.append(response[0].split('\">')[1])
        coupons.append(response[1])
        maturities.append(response[2])
        calldates.append(response[3])
        bids.append(response[4])
        asks.append(response[5])
        eods.append(response[6].split('</td>')[0])

    df = pd.DataFrame(
        {'CUSIP': cusips, 'SECURITY': types, "COUPON": coupons, "MATURITY": maturities, "CALLDATE": calldates,
         "BID": bids, "ASK": asks, "EOD": eods})

    return df


def query_type(cusip, price):
    with urllib.request.urlopen(
            "https://www.treasurydirect.gov/TA_WS/securities/search?cusip={}&format=json".format(cusip)) as url:
        data = url.read().decode()
        queries = dict()
        query = json.loads(data)
        queries[cusip] = query
    root = ET.Element('bond_data')
    node_cusip = ET.SubElement(root, 'cusip')
    node_cusip.text = cusip
    node_issue = ET.SubElement(root, 'issue')
    ET.SubElement(node_issue, "issuer").text = "US Treasury"
    ET.SubElement(node_issue, "issue_date").text = queries[cusip][0]['issueDate'][0:10]
    node_schedule = ET.SubElement(root, 'schedule')
    if queries[cusip][0]['interestRate'] == '':
        schedule = ET.SubElement(root, 'cusip')
    tree = ET.ElementTree(root)
    tree.write(os.path.join('data', '{}.xml'.format(cusip)))
    return queries[cusip], price


def bond_curve_fitter(instruments, calendar, daycounter=ql.ActualActual(), bdc=ql.ModifiedFollowing,
                      curve_method='LogCubicDiscount'):  # function to calibrate a spot curve
    bondSettlementDays, notional, redemption = 0, 100.0, 100.0
    bondHelper, depoHelper = [], []
    for tenor, instrument in instruments.items():
        maturity = tenor.split('-')
        reverse_maturity = ql.Period(-int(maturity[0]), units[maturity[1]])
        settlementDate = calendar.advance(instrument[2], reverse_maturity)
        schedule = ql.Schedule(settlementDate, instrument[2], freq[instrument[3]], calendar, bdc, bdc,
                               ql.DateGeneration.Backward, True)
        quote = ql.QuoteHandle(ql.SimpleQuote(float(instrument[5])))
        coupon = float(instrument[4]) if len(instrument[4]) > 0 else 0.0
        bondHelper.append(
            ql.FixedRateBondHelper(quote, bondSettlementDays, notional, schedule, [coupon / 100.0], daycounter, bdc))

    if curve_method == 'LogCubicDiscount':
        curve = ql.PiecewiseLogCubicDiscount(ql.Date.todaysDate(), bondHelper, daycounter)
    elif curve_method == 'LinearZero':
        curve = ql.PiecewiseLinearZero(ql.Date.todaysDate(), bondHelper, daycounter)
    elif curve_method == 'LinearFwd':
        curve = ql.PiecewiseLinearForward(ql.Date.todaysDate(), bondHelper, daycounter)
    elif curve_method == 'CubicZero':
        curve = ql.PiecewiseCubicZero(ql.Date.todaysDate(), bondHelper, daycounter)
    else:
        curve = ql.PiecewiseFlatForward(ql.Date.todaysDate(), bondHelper, daycounter)
        print('Curve Method unsupported')
    return curve, curve_method


def spot_curve_to_par_rate(curve, quote_dates=None, calc_date=ql.Date.todaysDate(),
                           daycounter=ql.ActualActual()):  # par rate approximation
    ts = ql.RelinkableYieldTermStructureHandle(curve)
    sum = 0.0
    dates = curve.dates() if quote_dates == None else quote_dates
    for i in range(1, len(dates)):
        yrs = daycounter.yearFraction(dates[i - 1], dates[i])
        sum += ts.discount(dates[i]) * yrs
    result = ts.discount(dates[0]) - ts.discount(dates[-1])
    return 100 * (result / sum)


def get_rates(curve, calc_date=ql.Date.todaysDate(), daycounter=ql.ActualActual()):
    spots, tenors, quote_dates = list(), list(), list()
    for d in curve.dates():
        yrs = daycounter.yearFraction(calc_date, d)
        compounding = ql.Compounded
        freq = ql.Semiannual
        zero_rate = curve.zeroRate(yrs, compounding, freq)
        quote_dates.append(d)
        tenors.append(yrs)
        eq_rate = zero_rate.equivalentRate(daycounter, compounding, freq, calc_date, d).rate()
        spots.append(100 * eq_rate)
    return quote_dates, tenors, spots


def bond_futures_analytics_with_CTD(curve, delivery_date, futPrice, bond_basket, calc_date=ql.Date.todaysDate(),
                                    calendar=ql.UnitedStates(), daycounter=ql.ActualActual()):
    cleanPrice = futPrice * curve.discount(delivery_date)
    settlement_days, face_value = 1, 100
    bondEngine = ql.DiscountingBondEngine(ql.YieldTermStructureHandle(curve))
    deliverables = list()
    min_basis, min_basis_index, cnt = 100., -1, 0
    for cusip, price in bond_basket.items():
        info = query_type(cusip, price)[0][0]
        issueDate, maturity, coupon = info['issueDate'][:-9], info['maturityDate'][:-9], float(
            info['interestRate']) / 100.0
        schedule = ql.Schedule(ql.DateParser.parseFormatted(issueDate, '%Y-%m-%d'),
                               ql.DateParser.parseFormatted(maturity, '%Y-%m-%d'), ql.Period(ql.Semiannual), calendar,
                               ql.ModifiedFollowing, ql.ModifiedFollowing, ql.DateGeneration.Forward, False)
        bond = ql.FixedRateBond(settlement_days, face_value, schedule, [coupon], daycounter)
        bond.setPricingEngine(bondEngine)
        zSpread = ql.BondFunctions.zSpread(bond, cleanPrice,
                                           curve, daycounter,
                                           ql.Compounded, ql.Semiannual,
                                           calc_date) * 10000
        conversionFactor = ql.BondFunctions.cleanPrice(bond, 0.06, daycounter, ql.Compounded, ql.Semiannual,
                                                       calc_date) / 100.
        adjFutPrice = futPrice * conversionFactor
        basis = price - adjFutPrice
        if basis < min_basis:
            min_basis, min_basis_index = basis, cnt
        print('ISIN: {}, zspread: {}bps'.format(cusip, zSpread))
        deliverables.append((cusip, bond, conversionFactor, price))
        cnt += 1
    print('Min Basis: {}'.format(min_basis))
    print('CTD Bond {} with conversion factor {} quoted at {}'.format(deliverables[min_basis_index][0],
                                                                      deliverables[min_basis_index][2],
                                                                      deliverables[min_basis_index][3]))
    futures = ql.FixedRateBondForward(calc_date, delivery_date, ql.Position.Long, 0.0, settlement_days, daycounter,
                                      calendar, ql.ModifiedFollowing, deliverables[min_basis_index][1],
                                      ql.YieldTermStructureHandle(curve), ql.YieldTermStructureHandle(curve))
    futuresTV = futures.cleanForwardPrice() / deliverables[min_basis_index][2]
    impliedYield = futures.impliedYield(deliverables[min_basis_index][3] / deliverables[min_basis_index][2], futPrice,
                                        calc_date, ql.Compounded, daycounter).rate()
    futZSpread = ql.BondFunctions.zSpread(deliverables[min_basis_index][1], deliverables[min_basis_index][3], curve,
                                          daycounter, ql.Compounded, ql.Semiannual, calc_date)
    fwdYield = ql.BondFunctions.bondYield(deliverables[min_basis_index][1], deliverables[min_basis_index][3],
                                          daycounter, ql.Compounded, ql.Semiannual, calc_date)

    print('Model Price {}, Market Price {}, Convexity Adjustment {}'.format(futuresTV, futPrice, futuresTV - futPrice))
    print('Implied Yield {}%'.format(impliedYield * 100))
    print('Forward ZSpread {}bps'.format(futZSpread * 10000))
    print('Forward YTM {}%'.format(fwdYield * 10000))


def credit_bond_price_to_zspread(bond, cleanprice, curve, method="Secant"):
    # implement bisection for a bond clean price to zspread
    tol, lower, upper = 1e-16, 0.0, 0.01
    ts = ql.RelinkableYieldTermStructureHandle(curve)
    engine = ql.DiscountingBondEngine(ts)
    bond.setPricingEngine(engine)
    assert (bond.NPV() > cleanprice)
    spreaded1 = ql.YieldTermStructureHandle(ql.ZeroSpreadedTermStructure(ts, ql.QuoteHandle(ql.SimpleQuote(upper))))
    engine1 = ql.DiscountingBondEngine(spreaded1)
    bond.setPricingEngine(engine1)
    if method == "Bisection":
        while (bond.NPV() > cleanprice):
            lower = upper
            upper = 2 * upper
            spreaded1 = ql.YieldTermStructureHandle(
                ql.ZeroSpreadedTermStructure(ts, ql.QuoteHandle(ql.SimpleQuote(upper))))
            engine1 = ql.DiscountingBondEngine(spreaded1)
            bond.setPricingEngine(engine1)
        print('Entering Bisection ...')
        mid = (lower + upper) / 2.0
        while (upper - lower) / 2.0 > tol:
            spreaded1 = ql.YieldTermStructureHandle(
                ql.ZeroSpreadedTermStructure(ts, ql.QuoteHandle(ql.SimpleQuote(lower))))
            engine1 = ql.DiscountingBondEngine(spreaded1)
            bond.setPricingEngine(engine1)
            fa = bond.NPV() - cleanprice
            spreaded1 = ql.YieldTermStructureHandle(
                ql.ZeroSpreadedTermStructure(ts, ql.QuoteHandle(ql.SimpleQuote(mid))))
            engine1 = ql.DiscountingBondEngine(spreaded1)
            bond.setPricingEngine(engine1)
            fc = bond.NPV() - cleanprice
            if abs(fc) < tol:
                return mid
            elif fa * fc < 0:
                upper = mid
            else:
                lower = mid
            mid = (lower + upper) / 2.0
        spreaded1 = ql.YieldTermStructureHandle(ql.ZeroSpreadedTermStructure(ts, ql.QuoteHandle(ql.SimpleQuote(upper))))
        engine1 = ql.DiscountingBondEngine(spreaded1)
        bond.setPricingEngine(engine1)
        # print('Price {}'.format(bond.NPV()))
        zspread = mid
    elif method == "Secant":
        guess0, guess1 = 0.0, 0.1
        bond.setPricingEngine(ql.DiscountingBondEngine(ql.YieldTermStructureHandle(
            ql.ZeroSpreadedTermStructure(ts, ql.QuoteHandle(ql.SimpleQuote(guess0))))))
        pv0 = bond.NPV()
        guess, pv = guess0, pv0
        print('Entering Quasi-Newton...')
        while (abs(guess0 - guess1) > tol):
            guess0, pv0 = guess, pv
            bond.setPricingEngine(ql.DiscountingBondEngine(ql.YieldTermStructureHandle(
                ql.ZeroSpreadedTermStructure(ts, ql.QuoteHandle(ql.SimpleQuote(guess1))))))
            pv1 = bond.NPV()
            if (abs(pv0 - pv1) < tol):
                break
            guess = guess1
            guess1 = guess1 - (pv1 - cleanprice) * (guess1 - guess0) / (pv1 - pv0)
            pv = pv1
        zspread = guess1
    else:
        zspread = None
    return zspread


def credit_bond_risk_measures(bond, cleanprice, curve):
    risks = list()
    bump_size = 1e-4
    bond_yield = bond.bondYield(cleanprice, ql.ActualActual(), ql.Compounded, ql.Semiannual)
    risks.append(('YieldToMaturity', bond_yield))
    rate = ql.InterestRate(bond_yield, ql.ActualActual(), ql.Compounded, ql.Semiannual)
    risks.append(('MacaulayDuration', ql.BondFunctions.duration(bond, rate, ql.Duration.Macaulay)))
    risks.append(('ModifiedDuration', ql.BondFunctions.duration(bond, rate, ql.Duration.Modified)))
    risks.append(('Convexity', ql.BondFunctions.convexity(bond, rate)))
    nodes = [ql.Period(1, ql.Months), ql.Period(3, ql.Months), ql.Period(6, ql.Months), ql.Period(1, ql.Years),
             ql.Period(2, ql.Years), ql.Period(3, ql.Years), ql.Period(5, ql.Years), ql.Period(10, ql.Years),
             ql.Period(15, ql.Years), ql.Period(30, ql.Years)]
    dates = [ql.Date.todaysDate() + node for node in nodes]
    buckets = ['1M', '3M', '6M', '1Y', '2Y', '3Y', '5Y', '10Y', '15Y', '30Y']
    discount_handle = ql.RelinkableYieldTermStructureHandle(curve)
    bond.setPricingEngine(ql.DiscountingBondEngine(discount_handle))
    base_price = bond.NPV()
    spreads = [ql.SimpleQuote(0.0) for i in range(len(nodes))]
    cnt = 0
    for spread in spreads:
        spread.setValue(bump_size)
        bumped_curve = ql.SpreadedLinearZeroInterpolatedTermStructure(ql.YieldTermStructureHandle(curve),
                                                                      [ql.QuoteHandle(q) for q in spreads], dates)
        discount_handle.linkTo(bumped_curve)
        bumped_price = bond.NPV()
        risks.append(('KeyRateDuration_{}'.format(buckets[cnt]), (bumped_price - base_price) / bump_size))
        spread.setValue(0.0)
        cnt += 1
    return risks


def callablebond_price_with_interest_lattice(callable, price, curve, reversion, volatility, grid_points):
    model = ql.HullWhite(ql.YieldTermStructureHandle(curve), reversion, volatility)
    engine = ql.TreeCallableFixedRateBondEngine(model, grid_points)
    callable.setPricingEngine(engine)
    cleanPrice = callable.cleanPrice()
    oas = callable.OAS(price, ql.YieldTermStructureHandle(curve), ql.ActualActual(), ql.Compounded, ql.Semiannual)
    return cleanPrice, oas


def test1():
    tenors = ['4-Week', '8-Week', '13-Week', '26-Week', '52-Week', '2-Year', '3-Year', '5-Year', '10-Year', '30-Year']
    calendar = ql.UnitedStates()
    quote_date = calendar.advance(ql.Date.todaysDate(), ql.Period(-1, ql.Days))  # get last business day EOD
    instruments = dict()
    for tenor in tenors:
        instruments[tenor] = ['', calendar.advance(quote_date,
                                                   ql.Period(-int(tenor.split('-')[0]), units[tenor.split('-')[1]])),
                              0.0]  # define earliest date for each calibrated instrument

    quote_date = dt.datetime(year=quote_date.year(), month=quote_date.month(), day=quote_date.dayOfMonth())
    data = grab_bond_data(quote_date)
    fittable = ['MARKET BASED BILL', 'MARKET BASED NOTE', 'MARKET BASED BOND']

    print('Test 0: Bootstrapping a daily yield curve using most recently auctioned T-instruments')
    print('Forming calibrated instruments...')
    pool = ThreadPoolExecutor(50)
    tasks = list()

    t1 = time.time()
    for (i, row) in data.iterrows():  # update latest auctioned instrument
        if row['SECURITY'] in fittable:
            tasks.append(pool.submit(query_type, row['CUSIP'], row['BID']))

    for t in tasks:
        quotes, price = t.result()
        for quote in quotes:
            if quote['securityTerm'] in instruments and dt.datetime.strptime(quote['auctionDate'][:-9],
                                                                             '%Y-%m-%d') > dt.datetime(
                year=instruments[quote['securityTerm']][1].year(),
                month=instruments[quote['securityTerm']][1].month(),
                day=instruments[quote['securityTerm']][1].dayOfMonth()):
                instruments[quote['securityTerm']] = [quote['cusip'],
                                                      ql.DateParser.parseFormatted(quote['auctionDate'][:-9],
                                                                                   '%Y-%m-%d)'),
                                                      ql.DateParser.parseFormatted(quote['maturityDate'][:-9],
                                                                                   '%Y-%m-%d)'),
                                                      quote['interestPaymentFrequency'], quote['interestRate'],
                                                      price]
    curve, curve_method = bond_curve_fitter(instruments, calendar)

    print(r'USD Daily Treasury curve completed in {} seconds'.format(time.time() - t1))
    observables, tenor, spot = get_rates(curve)
    spots, pars, cnt = list(), list(), 1
    for date in observables:
        if date > observables[0]:
            dates = observables[0:cnt + 1]
            spots.append((date, spot[cnt]))
            pars.append((date, spot_curve_to_par_rate(curve, dates)))
            cnt += 1
    print('Spot rates are {}'.format(spots))
    print('Par rates are {}'.format(pars))

    print('Test 1: Calculating measures for Bond Futures')
    print('Forming deliverable bonds')
    bondBasket = {'912828U32': 99.140625, '9128283H1': 99.5625, '912828G61': 99.40625, '912828UB4': 99.09375,
                  '912828U73': 99.28125, '9128283N8': 99.59375, '912828UF5': 99.0625, '912828G95': 99.40625}
    delivery_date, futPrice = ql.IMM_date('M9'), 115 + 13.25 / 32
    bond_futures_analytics_with_CTD(curve, delivery_date, futPrice, bondBasket)

    print('\nTest2: Calculating credit spread and risk measures for a corporate bond')
    isin, price = 'US852061AA81', 99.00
    print('ISIN: {} quoted at {}'.format(isin, price))
    coupon, issue_date, maturity = 0.0925, ql.DateParser.parseFormatted('1992-04-15',
                                                                        '%Y-%m-%d)'), ql.DateParser.parseFormatted(
        '2022-04-15', '%Y-%m-%d)')
    schedule = ql.Schedule(issue_date, maturity, ql.Period(ql.Semiannual), ql.UnitedStates(), ql.ModifiedFollowing,
                           ql.ModifiedFollowing, ql.DateGeneration.Backward, False)
    bond = ql.FixedRateBond(1, 100.0, schedule, [coupon], ql.ActualActual(), ql.ModifiedFollowing)
    t1 = time.time()
    print('Credit spread is {} for quoted price {}, completed in {} seconds'.format(
        credit_bond_price_to_zspread(bond, price, curve, 'Bisection'), price, time.time() - t1))
    t1 = time.time()
    print('Credit spread is {} for quoted price {}, completed in {} seconds'.format(
        credit_bond_price_to_zspread(bond, price, curve, 'Secant'), price, time.time() - t1))

    print('Calculating risk measures for a corporate bond')
    t1 = time.time()
    risks = credit_bond_risk_measures(bond, price, curve)
    for val in risks:
        print('Risk measure {} has value {}'.format(val[0], val[1]))
    print('Risk measures completed in {} seconds'.format(time.time() - t1))

    print('\nTest 3: Calculating option adjusted spread and price with short-rate model for a callable bond')
    isin, price = 'US06060WBJ36', 97.15
    print('ISIN: {} quoted at {}'.format(isin, price))
    coupon, issue_date, maturity = 0.0465, ql.DateParser.parseFormatted('2012-04-15',
                                                                        '%Y-%m-%d)'), ql.DateParser.parseFormatted(
        '2022-04-15', '%Y-%m-%d)')
    schedule = ql.Schedule(issue_date, maturity, ql.Period(ql.Quarterly), ql.UnitedStates(), ql.ModifiedFollowing,
                           ql.ModifiedFollowing, ql.DateGeneration.Backward, False)
    bond = ql.FixedRateBond(1, 100.0, schedule, [coupon], ql.ActualActual(), ql.ModifiedFollowing)

    callschedule = ql.CallabilitySchedule()
    for cashflow in bond.cashflows():
        if cashflow.date() > ql.Date(15, 4, 2020):
            callschedule.append(
                ql.Callability(ql.CallabilityPrice(100.0, ql.CallabilityPrice.Clean), ql.Callability.Call,
                               cashflow.date()))
    callable = ql.CallableFixedRateBond(1, 100.0, schedule, [coupon], ql.ActualActual(), ql.ModifiedFollowing, 100.0,
                                        issue_date, callschedule)
    reversion, volatility, grid = 0.03, 0.01, 40
    print('Build a Hull-White with reversion {} and vol {}, preset grid step {}'.format(reversion, volatility, grid))
    t1 = time.time()
    cleanPrice, optionAdjustedSpread = callablebond_price_with_interest_lattice(callable, price, curve, reversion,
                                                                                volatility, grid)
    print(
        "Callable clean price is {} with short-rate model setting, OAS is {} for quoted price {}, completed in {} seconds".format(
            cleanPrice, optionAdjustedSpread, price, time.time() - t1))
    reversion, volatility, grid = 0.03, 0.01, 100
    print('Build a Hull-White with reversion {} and vol {}, preset grid step {}'.format(reversion, volatility, grid))
    t1 = time.time()
    cleanPrice, optionAdjustedSpread = callablebond_price_with_interest_lattice(callable, price, curve, reversion,
                                                                                volatility, grid)
    print(
        "Callable clean price is {} with short-rate model setting, OAS is {} for quoted price {}, completed in {} seconds".format(
            cleanPrice, optionAdjustedSpread, price, time.time() - t1))


if __name__ == '__main__':
    test1()