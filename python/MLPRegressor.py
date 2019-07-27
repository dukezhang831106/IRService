import pandas as pd
import numpy as np
import keras
from keras.models import Sequential, Model
from keras.layers import Dense, Input
from keras import optimizers
from keras.optimizers import Adam

import seaborn as sns
import matplotlib.pyplot as plt
sns.set()


class CurveDynamics(object):
    def __init__(self, data):
        sns.set(style="ticks")
        self.data = data
        self.tenors = ['1Y', '2Y', '3Y', '4Y', '5Y', '6Y', '7Y', '8Y', '9Y', '10Y', '15Y', '20Y', '30Y']
        self.t = [float(tenor[:-1]) for tenor in self.tenors]
        self.data = self.data.set_index('DATE').dropna()
        self.abs_return = data.set_index('DATE').diff().dropna()

    def plot_raw_data(self):
        fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(15, 15))
        # plot 1: swap rate with tenor vs. trade date
        for tenor in self.tenors:
            axes[0, 0].plot(self.data.index.values, self.data[tenor])
        axes[0, 0].set_xlabel('Date')
        axes[0, 0].set_ylabel('Rate')
        axes[0, 0].set_xticks(range(0, len(self.data), 200), minor=False)
        axes[0, 0].set_title('Swap Rate')

        # plot 2: box and whiskier: distribution vs. tenor
        dist = self.data.transpose()
        axes[0, 1].boxplot(dist, whis=1, positions=[int(tenor[:-1]) for tenor in self.tenors], showfliers=False)
        axes[0, 1].set_xlabel('Tenor')
        axes[0, 1].set_ylabel('Rate')
        axes[0, 1].grid(True)
        axes[0, 1].set_title('Swap Rate Distribution')

        # plot 3: absolute return
        for tenor in self.tenors:
            axes[1, 0].plot(self.abs_return.index.values, self.abs_return[tenor])
        axes[1, 0].set_xlabel('Date')
        axes[1, 0].set_ylabel('Return')
        axes[1, 0].set_xticks(range(0, len(self.data), 200), minor=False)
        axes[1, 0].set_title('Swap Rate Absolute Return')


        # plot 4: return correlation
        corr = self.abs_return[self.tenors].corr()
        sns.heatmap(corr, ax=axes[1,1])

        plt.suptitle('Swap Rate statistics', fontsize=18)
        plt.show()

    def PCARegressor(self, n=3, train_ratio = 0.8):
        from sklearn.decomposition import PCA
        import numpy as np
        train_size = int(len(self.data) * train_ratio)

        X = np.array(self.data.head(train_size)[self.tenors])
        self.pca = PCA(n_components=n)
        self.pca.fit(X)

        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 5))
        axes[0].bar(range(n), self.pca.explained_variance_ratio_)
        plt.sca(axes[0])
        for i, v in enumerate(self.pca.explained_variance_ratio_):
            axes[0].text(i-0.125, v, '{}%'.format(round(100*v,2)))
        plt.xticks(range(n), ['Component {}'.format(i+1) for i in range(n)])
        vals = axes[0].get_yticks()
        axes[0].set_yticklabels(['{:,.2%}'.format(x) for x in vals])
        axes[0].set_xlabel('Components')
        axes[0].set_ylabel('Variance Ratio')
        axes[0].set_title('Variance explained by Principal Components')

        for component in self.pca.components_:
            axes[1].plot(self.t, component)
        axes[1].set_xlabel('Tenor')
        axes[1].set_ylabel('Weight')
        axes[1].set_title('Principal Components Eigenvectors')

        plt.suptitle('PCA Analysis', fontsize=18)
        plt.show()

        return

    def PCAShock(self, pivot = '5Y', cnt = 1, shock_sizes = [-1, 1]):
        from sklearn.linear_model import LinearRegression
        data = self.data
        data['Component_0'] = self.abs_return[pivot]
        components = np.matrix(data[self.tenors])*self.pca.components_.transpose()
        for i in range(self.pca.n_components_):
            data['Component_{}'.format(i+1)] = components[:, i]
        tags = ['Component_{}'.format(i) for i in range(self.pca.n_components_ + 1)]
        data = data[tags].join(self.abs_return, how='outer').dropna()
        model = LinearRegression()
        X, Y = np.array(data[tags]), np.array(data[self.tenors])
        model.fit(X, Y)

        base = np.matrix(self.data.iloc[cnt][self.tenors])
        df = pd.DataFrame({'Tenor': self.tenors, 'Base': [base[0, i] for i in range(len(self.tenors))]})
        plt.plot('Tenor', 'Base', data=df)
        for shock_size in shock_sizes:
            shockX = np.matrix(data.iloc[cnt][tags])
            shockX[0, 0] = shock_size
            shockY = model.predict(shockX)
            shocked = base + shockY
            df['PCA_{}bp'.format(100*shock_size)] = [shocked[0, i] for i in range(len(self.tenors))]
            plt.plot('Tenor', 'PCA_{}bp'.format(100*shock_size), data=df, marker='*')

        plt.legend()
        plt.xlabel('Tenor')
        plt.ylabel('Rate')
        plt.title('PCA Shock')
        plt.show()



    def MLP(self, pivot = '5Y', layers = 2, train_ratio = 0.8):
        num_tenors = len(self.data.transpose())
        input_size = num_tenors + 1
        hidden_size = num_tenors + 2
        output_size = num_tenors
        train_size = int(len(self.data) * train_ratio)

        data = Input(shape=(input_size,))
        inflow = data

        for layer in range(layers):
            outflow = Dense(hidden_size, activation="tanh")(inflow)
            inflow = outflow

        outflow = Dense(output_size, activation="tanh")(inflow)
        output = Dense(output_size, activation="linear")(outflow)

        model = Model(input=data, output=output)
        model.compile(optimizer='adam', loss='mse')
        data = self.data
        data['Dates'] = data.index.values
        ddata = self.abs_return
        ddata['Dates'] = ddata.index.values
        train_set = pd.merge(data, ddata, on='Dates', suffixes=('_Rate', '_Delta')).set_index('Dates')
        Xlabel = ['{}_Rate'.format(tenor) for tenor in self.tenors]
        Xlabel.append('{}_Delta'.format(pivot))
        Ylabel = ['{}_Delta'.format(tenor) for tenor in self.tenors]
        train_set_X = np.array(train_set[Xlabel].head(train_size))
        train_set_Y = np.array(train_set[Ylabel].head(train_size))
        train_set_X = train_set_X.reshape((train_size, input_size))
        train_set_Y = train_set_Y.reshape((train_size, output_size))

        model.fit(train_set_X, train_set_Y, epochs=3200, verbose=0)

        test_set_X = np.array(train_set[Xlabel].tail(len(self.data) - train_size))
        test_set_Y = np.array(train_set[Ylabel].tail(len(self.data) - train_size))
        test_set_Yhat = model.predict(test_set_X, verbose=0)


        for i in range(1):
            plt.plot(self.t, test_set_Y[i])
            plt.plot(self.t, test_set_Yhat[i], '+')
        plt.show()

        return test_set_Y, test_set_Yhat

    def AutoEncoder(self, layers = 1, train_ratio = 0.8):
        hidden_layer = layers
        input_size = len(self.data.transpose())
        hidden_size = int(input_size/4)
        output_size = input_size
        train_size = int(len(self.data)*train_ratio)

        data = Input(shape=(input_size, ))
        encoder = Dense(hidden_size, activation="tanh")(data)
        decoder = Dense(output_size, activation="tanh")(encoder)
        output = Dense(output_size, activation="linear")(decoder)

        self.autoencoder = Model(input=data, output=output)
        self.autoencoder.compile(optimizer='adam', loss='mse')
        X = np.array(self.data.head(train_size))
        X = X.reshape((train_size, input_size))
        self.autoencoder.fit(X, X, epochs=1600, verbose=0)

        return

    def lowDimensionRepresentation(self, cnt = 0):
        Data = np.array(self.data)
        DataHat = self.autoencoder.predict(Data)
        df = pd.DataFrame({'Tenor': self.t, 'exact': Data[cnt], 'lowDim': DataHat[cnt]})
        plt.plot('Tenor', 'exact', data = df)
        plt.plot('Tenor', 'lowDim', data = df, marker='*')
        plt.legend()
        plt.xlabel('Tenor')
        plt.ylabel('Rate')
        plt.title('Exact Data vs. AutoEncoder')
        plt.show()

        return

if __name__ == '__main__':
    data = pd.read_csv('/home/e658880/USDSwap14.csv')
    curve = CurveDynamics(data)
    curve.plot_raw_data()
    curve.PCARegressor()
    curve.AutoEncoder()
    curve.lowDimensionRepresentation(100)
    curve.PCAShock()
    #print(X)
