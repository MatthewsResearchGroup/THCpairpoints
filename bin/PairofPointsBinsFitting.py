# import related libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
import matplotlib.ticker as mtick
import os


path = os.getcwd()
csvpath = path + "/bootstrapping.csv"

figurepath = path + "/binfitting/"
if (not os.path.exists(figurepath)):
    os.makedirs(figurepath)

# load the data
df=pd.read_csv(csvpath)


# setting type

moles=df["mol"].unique()
methods=df["method"].unique()
basises=df["basis"].unique()
Ftypes=df["Ftype"].unique()
energy_types = df["energy_type"].unique()
ratios = df["# * (nv+no)"].unique()


def powerorder1(x, a ):
    return np.log10(a * a / x )


def powerorder2(x, a):
    return np.log10(a * a / (x * x) )

def powerorder3(x, a):
    return np.log10(a * a / (x * x * x) )

def powerorder1and2(x, a, b):
    return np.log10(a * a / x + b * b / (x*x))

def biexpo(x, a, b, c, d):
    return np.log10(a * a * np.exp(b * x) + c * c * np.exp(d * x))


def expo(x, a, b):
    return np.log10(a * a * np.exp(b * x))


#para_list  = []
for mol in moles:
    for method in methods:
        for basis in basises:
            fig, ax = plt.subplots(figsize=(10, 6), dpi = 600)
            for i in range(50):
                # here we should use dataframe
                savefigurepath =  figurepath + mol + '_' + method + '_' + basis + '.jpg'
                Ftype = 'all'
                y_name = "ave" + str(i+1)
                bin_name = "bin" + str(i+1)
                ex_type = 'Ex_add/Ex+Ec'
                ec_type = 'Ec_add/Ex+Ec'
                data_ex = df[(df['energy_type'] == ex_type) & (df['method'] == method) &(df['basis']== basis) & (df['Ftype'] == Ftype)][["# * (nv+no)", y_name, bin_name]]
                data_ec = df[(df['energy_type'] == ec_type) & (df['method'] == method) &(df['basis']== basis) & (df['Ftype'] == Ftype)][["# * (nv+no)", y_name, bin_name]]
            #     data_x_y = data_ex[["# * (nv+no)", y_name, bin_name]]
                data_x_y = data_ex[data_ex[y_name] < 0.0]
                if (data_x_y[y_name].sum() >= 0.0):
                    continue
                y_data = np.log10(-np.array(data_x_y[y_name]))
                x = np.array(data_x_y["# * (nv+no)"])
                distance = np.array(data_x_y[bin_name])[0]
                popt_x_bi, pcov_x_bi = curve_fit(biexpo, x, y_data ,p0=(1.0, 1.0, 1.0, 1.0), maxfev=5000)
                y_data_bi = biexpo(x, *popt_x_bi)
                R2_bi = r2_score(y_data, y_data_bi)
                #if (abs(popt_x_bi[1]) > abs(popt_x_bi[3])):
                #    para_list.append([distance, popt_x_bi[0] * popt_x_bi[0], popt_x_bi[1], popt_x_bi[2] * popt_x_bi[2], popt_x_bi[3], R2_bi])
                #else:
                #    para_list.append([distance, popt_x_bi[2] * popt_x_bi[2], popt_x_bi[3], popt_x_bi[0] * popt_x_bi[0], popt_x_bi[1], R2_bi])
                if (i < 6):
                    ax.scatter(x, y_data, label = "x", s=3)
                    # plt.scatter(x, y_c, label = "c", s=3, marker = '^')
                    # plt.plot(x, y_data_power1 , '-', label = 'power1, R2 = {}'.format(R2_power1))
                    # plt.plot(x, y_data_power2 , '-', label = 'power2, R2 = {}'.format(R2_power2))
                    # plt.plot(x, y_data_power3 , '-', label = 'power3, R2 = {}'.format(R2_power3))
            #         plt.plot(x, y_data_power12 , '-', label = 'power12, R2 = {}'.format(R2_power12))
            #         plt.plot(x, y_data_s, '-', label = 'expo, R2 = {}'.format(R2_s))
                    ax.plot(x, y_data_bi , '-', label = 'biexpo, R2 = {}'.format(R2_bi))
            ax.set_title('{} the first 6th bin, {} {}'.format(mol, method, basis))
            #         plt.legend()
            #         plt.show()a
            plt.savefig(savefigurepath, dpi = 600)
