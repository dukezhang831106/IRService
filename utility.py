import eikon as ek
import itertools
import datetime as dt
ek.set_app_key('3edad9e062734052908433f7c74ede994f827026')

def generate_queries(quote_date, ccy, instruments, instrument, sources):
    import QuantLib as ql
    commands = list()
    source = sources[instrument]
    val = instruments[instrument]
    if instrument == 'Depo':
        expiries = ['TN', 'SN', '1M', '3M']
        instruments = ['USD{}{}='.format(expiry, val) for expiry in expiries]
    elif instrument == 'Fut':
        expiries = list()
        current_date = ql.DateParser_parse(quote_date, '%Y-%m-%d')
        for i in range(8):
            expiries.append(ql.IMM_nextCode(current_date))
            current_date = ql.IMM_nextDate(current_date)
        instruments = ['{}{}'.format(val, expiry) for expiry in expiries]
    elif instrument == 'Swap':
        expiries = ['3Y', '4Y', '5Y', '6Y', '7Y', '8Y', '9Y', '10Y', '15Y', '20Y', '25Y', '30Y']
        instruments = ['USD{}{}={}'.format(val, expiry, source) for expiry in expiries]
    elif instrument == 'Swaption':
        expiries = ['1M', '3M', '6M', '1Y', '2Y', '3Y', '4Y', '5Y', '6Y', '7Y', '8Y', '9Y', '10Y', '15Y', '20Y', '25Y', '30Y']
        tenors = ['1Y', '2Y', '3Y', '4Y', '5Y', '6Y', '7Y', '8Y', '9Y', '10Y', '15Y', '20Y', '25Y', '30Y']
        instruments = ['USD{}{}X{}={}'.format(val, expiry, tenor, source) for expiry in expiries for tenor in tenors]

    result, err = ek.get_data(instruments, ['CF_LAST'])

    #instruments = ['USD1YX1YN{}=R'.format(i+1) for i in range(13)]
    #result, err = ek.get_data(instruments, ['CF_LAST'])
    #print(result)
    for i in range(len(result)):
        if instrument != 'Swaption':
            expiry, tenor, quote = expiries[i], '', result.iloc[i, :]['CF_LAST']
        else:
            bucket = result.iloc[i, :]['Instrument'].split('=')[0][3:].split('X')
            expiry, tenor, quote = bucket[0], bucket[1], result.iloc[i, :]['CF_LAST']
        commands.append("INSERT INTO \"IRInstruments\" VALUES (Date('{}'), '{}', '{}', '{}', '{}', '{}', {})".format(quote_date, ccy, instrument, source, expiry, tenor, quote))
    return commands

def update_database(commands):
    import psycopg2
    from psycopg2.extensions import ISOLATION_LEVEL_AUTOCOMMIT
    connection = psycopg2.connect(dbname='EikonInstrument', user='postgres', host='127.0.0.1', port='1111', password='123456')
    cursor = connection.cursor()
    #cursor.execute("ALTER TABLE \"IRInstruments\" ADD COLUMN Currency TEXT")
    #cursor.execute("ALTER TABLE \"IRInstruments\" ADD COLUMN InstrumentType TEXT ")
    #cursor.execute("ALTER TABLE \"IRInstruments\" ADD COLUMN Source TEXT ")
    #cursor.execute("ALTER TABLE \"IRInstruments\" ADD COLUMN Expiry TEXT ")
    #cursor.execute("ALTER TABLE \"IRInstruments\" ADD COLUMN Tenor TEXT ")
    #cursor.execute("ALTER TABLE \"IRInstruments\" ADD COLUMN Quote double precision ")

    for command in commands:
        cursor.execute(command)
    connection.commit()
    return

def parse_database_to_csv(path):
    import psycopg2
    import pandas as pd
    connection = psycopg2.connect(dbname='EikonInstrument', user='postgres', host='127.0.0.1', port='1111', password='123456')
    cursor = connection.cursor()
    cursor.execute("SELECT * FROM \"IRInstruments\"")
    instruments = cursor.fetchall()
    dates, ccys, types, sources, expiries, tenors, quotes = list(), list(), list(), list(), list(), list(), list()
    for instrument in instruments:
        dates.append(instrument[0])
        ccys.append(instrument[1])
        types.append(instrument[2])
        sources.append(instrument[3])
        expiries.append(instrument[4])
        tenors.append(instrument[5])
        quotes.append(instrument[6])
    df = pd.DataFrame({'Dates': dates, 'Currency': ccys, 'Instruments': types, 'Sources': sources, 'Expiry': expiries, 'Tenor': tenors, 'Quote': quotes})
    df.to_csv(path, index=False)
    connection.close()



if __name__ == "__main__":
    #update_database('')
    instruments = {'Depo':'D', 'Fut':'ED', 'Swap':'AM3L', 'Swaption':''}
    quote_date = dt.date.today().strftime('%Y-%m-%d')
    ccys = ['USD']
    sources = {'Depo':'', 'Fut':'', 'Swap':'ICAP', 'Swaption':'ICAP'}
    #keys = ['Swaption']
    keys = ['Depo', 'Fut', 'Swap', 'Swaption']
    for ccy in ccys:
        for key in keys:
            commands = generate_queries(quote_date, ccy, instruments, key, sources)
            update_database(commands)
            print('{} in {} saved'.format(key, ccy))


#data, error = ek.get_data('.TRISTI', ['TR.CLOSEPRICE.Date', 'TR.CLOSEPRICE'], parameters={'SDate':'2019-01-01', 'EDate':'2019-05-30'})
#print(data)