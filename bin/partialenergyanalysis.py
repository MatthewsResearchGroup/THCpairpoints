# import related libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
import os
import math
np.seterr(divide = 'ignore') 

# read the data and get some basic information like molecule name, method, and basis set
path="./partialGridEnergy.csv"
df=pd.read_csv(path)

molecules=df["Molecule"].unique()
methods=df["Method"].unique()
basises=df["Basis"].unique()

# add more analysis.
df["ratio_grid_num"]= df["num_grid_keep"] / df["norb"]
df["PerError_c"] = abs( (df["Ec_THC"] - df["Ec_exact"]) / df["Ec_exact"])
df["PerError_x"] = abs( (df["Ex_THC"] - df["Ex_exact"]) / df["Ex_exact"])
df["PerError_xminsc"] = df["PerError_x"] - df["PerError_c"]

# get rid of the non-positive value since we need to take the log of PerError_x or PerError_c
df  = df[(df['PerError_c'] > 0) & (df['PerError_x'] > 0)]


# define some fitting functions here, you can add your own fitting functions.
def biexpo(x, a, b, c, d):
    return np.log10(a * a * np.exp(b * x) + c * c * np.exp(d * x))

def biexpofixedexpo(x, a, c):
    global b, d
    return np.log10(biexpo(x, a, bglobal, c, dgolbal))

def expo(x, a, b):
    return np.log10(a * a * np.exp(b * x))

def inversepower1(x, a):
    return np.log10(a * a / x)

def inversepower2(x, a):
    return np.log10(a * a / (x * x) )

def inversepower3(x, a):
    return np.log10(a * a / (x * x * x) )

def biinversepower(x, a, b, c, d):
    return np.log10(a*a/(x**b) + c*c/ (x**d))


d = ['Molecule', 'Method', 'Basis', 'grid_point_Ec_lt_0.1p', 'Ex_when_Ec_lt_0.1p', 'Ex_biexpo_a', 'Ex_biexpo_b', 'Ex_biexpo_c', 'Ex_biexpo_d', 'Ex_biexpo_R2', 'Ex_expo_a', 'Ex_expo_b','Ex_expo_R2', 'Ec_biexpo_a', 'Ec_biexpo_b', 'Ec_biexpo_c', 'Ec_biexpo_d', 'Ec_biexpo_R2', 'Ec_expo_a', 'Ec_expo_b','Ec_expo_R2']
df_save = pd.DataFrame(columns=d)


for mol in molecules:
    for method in methods:
        for basis in basises:
            print(mol, method, basis)
            data = df[(df['Molecule'] == mol) & (df['Method'] == method) & (df['Basis'] == basis)]
            x = data["ratio_grid_num"]
            y_x = np.log10(data["PerError_x"])
            y_c = np.log10(data["PerError_c"]) 
            # regression for exchange energy
            #             popt_x_bi, pcov_x_bi = curve_fit(biexpo,x,y_x,p0=(1.0,-1.0,1.0,-1.0), maxfev=5000)
            popt_x_s, pcov_x_s = curve_fit(expo,x,y_x,p0=(1.0,-1.0), maxfev=5000)
            popt_x_power1, pcov_x_power1 = curve_fit(inversepower1,x,y_x,p0=(0.1), maxfev=500)
            popt_x_power2, pcov_x_power2 = curve_fit(inversepower2,x,y_x,p0=(0.1), maxfev=500)
            popt_x_power3, pcov_x_power3 = curve_fit(inversepower3,x,y_x,p0=(0.1), maxfev=500)
            popt_x_bipower, pcov_x_bipower = curve_fit(biinversepower,x,y_x,p0=(0.1, 1, 0.1, 2), maxfev=500)


            # initial guess for the biexpo which use the first 20 and last 20 point to fit
            popt_x_bi_f20, pcov_x_bi_f20 = curve_fit(expo,x[0:30],y_x[0:30],p0=(1.0,-1.0), maxfev=5000)
            popt_x_bi_l20, pcov_x_bi_l20 = curve_fit(expo,x[-30:],y_x[-30:],p0=(1.0,-1.0), maxfev=5000)

            # fit the coeff part 


            popt_x_bi, pcov_x_bi = curve_fit(biexpo,x,y_x,p0=(1.0, popt_x_bi_f20[1],1.0,popt_x_bi_l20[1]), maxfev=5000)
            y_x_bi = biexpo(x, *popt_x_bi)
            y_x_s = expo(x, *popt_x_s)
            y_x_power1 = inversepower1(x, *popt_x_power1)
            y_x_power2 = inversepower2(x, *popt_x_power2)
            y_x_power3 = inversepower3(x, *popt_x_power3)
            y_x_bipower = biinversepower(x, *popt_x_bipower)

            #             y_x_s = expo(x, *popt_x_s)
            R2_x_bi = r2_score(y_x, y_x_bi)
            R2_x_s = r2_score(y_x, y_x_s)
            R2_x_power1 = r2_score(y_x, y_x_power1)
            R2_x_power2 = r2_score(y_x, y_x_power2)
            R2_x_power3 = r2_score(y_x, y_x_power3)
            R2_x_bipower = r2_score(y_x, y_x_bipower)

            # regression for columb energy
            #             popt_c_bi, pcov_c_bi = curve_fit(biexpo,x,y_c,p0=(1.0,-1.0,1.0,-1.0), maxfev=5000)
            popt_c_s, pcov_c_s = curve_fit(expo,x,y_c,p0=(1.0,-1.0), maxfev=5000)
            popt_c_power1, pcov_c_power1 = curve_fit(inversepower1,x,y_c,p0=(0.1), maxfev=500)
            popt_c_power2, pcov_c_power2 = curve_fit(inversepower2,x,y_c,p0=(0.1), maxfev=500)
            popt_c_power3, pcov_c_power3 = curve_fit(inversepower3,x,y_c,p0=(0.1), maxfev=500)
            popt_c_bipower, pcov_c_bipower = curve_fit(biinversepower,x,y_c,p0=(0.1, 1, 0.1, 2), maxfev=500)


            popt_c_bi_f20, pcov_c_bi_f20 = curve_fit(expo,x[0:30],y_c[0:30],p0=(1.0,-1.0), maxfev=5000)
            popt_c_bi_l20, pcov_c_bi_l20 = curve_fit(expo,x[-30:],y_c[-30:],p0=(1.0,-1.0), maxfev=5000)

            popt_c_bi, pcov_c_bi = curve_fit(biexpo,x,y_c,p0=(1.0,popt_c_bi_f20[1],1.0,popt_c_bi_l20[1]), maxfev=5000)

            # get fiited y
            y_c_bi = biexpo(x, *popt_c_bi)
            y_c_s = expo(x, *popt_c_s)
            y_c_power1 = inversepower1(x, *popt_c_power1)
            y_c_power2 = inversepower2(x, *popt_c_power2)
            y_c_power3 = inversepower3(x, *popt_c_power3)
            y_c_bipower = biinversepower(x, *popt_c_bipower)


            # calculate R2
            R2_c_bi = r2_score(y_c, y_c_bi)
            R2_c_s = r2_score(y_c, y_c_s) 
            R2_c_power1 = r2_score(y_c, y_c_power1)
            R2_c_power2 = r2_score(y_c, y_c_power2)
            R2_c_power3 = r2_score(y_c, y_c_power3) 
            R2_c_bipower = r2_score(y_c, y_c_bipower)

            turning_point_x = np.log(popt_x_bi[0]*popt_x_bi[0]/(popt_x_bi[2]*popt_x_bi[2]))/(popt_x_bi[3] - popt_x_bi[1])
            turning_point_c = np.log(popt_c_bi[0]*popt_c_bi[0]/(popt_c_bi[2]*popt_c_bi[2]))/(popt_c_bi[3] - popt_c_bi[1])
            errorx_tp_x = 10**(biexpo(turning_point_x, *popt_x_bi))
            errorx_tp_c = 10**(biexpo(turning_point_c, *popt_x_bi))
            errorc_tp_x = 10**(biexpo(turning_point_x, *popt_c_bi))
            errorc_tp_c = 10**(biexpo(turning_point_c, *popt_c_bi))

            # get the first point whose percentage error of coulomb energy is less than 0.1%, we need get how many 
            # grid points we use and the exchange energy percentage error when percentage error of coulomb energy is less than 0.1%
            data_c_less_than_1_over_1000 = data[data['PerError_c'] < 0.001]
            grid_c_less_than_1_over_1000 = data_c_less_than_1_over_1000['ratio_grid_num'].iloc[0]
            Ex_when_c_less_than_1_over_1000 = data_c_less_than_1_over_1000['PerError_x'].iloc[0]
            temp_list = [mol, method, basis,  grid_c_less_than_1_over_1000, Ex_when_c_less_than_1_over_1000, popt_x_bi[0]*popt_x_bi[0], popt_x_bi[1], popt_x_bi[2]*popt_x_bi[2], popt_x_bi[3],  R2_x_bi, popt_x_s[0]*popt_x_s[0], popt_x_s[1], R2_x_s, popt_c_bi[0]*popt_c_bi[0], popt_c_bi[1], popt_c_bi[2]*popt_c_bi[2], popt_c_bi[3],  R2_c_bi, popt_c_s[0]*popt_c_s[0], popt_c_s[1], R2_c_s]
            df_save.loc[len(df_save.index)] =  temp_list 

            # plotting the fitting figure. 

            path = "./partial_energy_fitting/" + mol
            isExist = os.path.exists(path)
            if (not isExist):
                os.makedirs(path) 

            figname = path + "/" + mol + "_" + method.strip() + "_" + basis + ".jpg"
            #             fig = plt.figure(figsize=(10, 6), dpi=600)
            #             figs, axs = plt.subplots(nrows=2, ncols=2)
            figs, axs = plt.subplots(nrows=1, ncols=2, figsize=(15, 5), dpi=600)
            plt.subplots_adjust(wspace=0.25)
            plt.rcParams['axes.linewidth'] = 2

            axs[1].tick_params(axis="y",direction="in", width = 1.5, length = 5, labelsize=15)
            axs[1].tick_params(axis="x",direction="in", width = 1.5, length = 5, labelsize=15)
            axs[0].tick_params(axis="y",direction="in", width = 1.5, length = 5, labelsize=15)
            axs[0].tick_params(axis="x",direction="in", width = 1.5, length = 5, labelsize=15)


            x_floor_plus1 = math.floor(max(x)) + 1
            x_new =  np.linspace(0, x_floor_plus1)
            #                 print("Here we go")
            axs[1].scatter(x, y_x, label = "x", s=4)
            axs[1].scatter(x, y_c, label = "c", s=4, marker = '^')
            axs[1].plot(x_new, biexpo(x_new, *popt_x_bi), '-', label='x_biexpo_fitting', lw = 2.5)
            axs[1].plot(x_new, expo(x_new, *popt_x_s), '-', label='x_expo_fitting', lw = 2.5)
            axs[1].plot(x_new, biexpo(x_new, *popt_c_bi), '-', label='c_biexpo_fitting', lw = 2.5)
            axs[1].plot(x_new, expo(x_new, *popt_c_s), '-', label='c_expo_fitting', lw = 2.5)

            # inverse power series fitting
            #             plt.plot(x_new, inversepower1(x_new, *popt_x_power1), '-', label='ex a/x_fitting: y = {:.4f} / x \n. R^2 = {:.4f}'.format(popt_x_power1[0] * popt_x_power1[0], R2_x_power1))
            #             plt.plot(x_new, inversepower2(x_new, *popt_x_power2), '-', label='ex a/x^2_fitting: y = {:.4f} / x^2 \n. R^2 = {:.4f}'.format(popt_x_power2[0] * popt_x_power2[0], R2_x_power2))
            #             plt.plot(x_new, inversepower3(x_new, *popt_x_power3), '-', label='ex a/x^3_fitting: y = {:.4f} / x^3 \n. R^2 = {:.4f}'.format(popt_x_power3[0] * popt_x_power3[0], R2_x_power3))
            #            plt.plot(x_new, biinversepower(x_new, *popt_x_bipower), '-', label='ex bi inverse power fitting: y = {:.4f} / x^{:.4f} + {:.4f} / x^{:.4f}  \n. R^2 = {:.4f}'.format(popt_x_bipower[0] * popt_x_bipower[0], popt_x_bipower[1], popt_x_bipower[2] * popt_x_bipower[2], popt_x_bipower[3] , R2_x_bipower))
            #             plt.plot(x_new, inversepower1(x_new, *popt_c_power1), '-', label='ec a/x_fitting: y = {:.4f} / x \n. R^2 = {:.4f}'.format(popt_c_power1[0] * popt_c_power1[0], R2_c_power1))
            #             plt.plot(x_new, inversepower2(x_new, *popt_c_power2), '-', label='ec a/x^2_fitting: y = {:.4f} / x^2 \n. R^2 = {:.4f}'.format(popt_c_power2[0] * popt_c_power2[0], R2_c_power2))
            #             plt.plot(x_new, inversepower3(x_new, *popt_c_power3), '-', label='ec a/x^3_fitting: y = {:.4f} / x^3 \n. R^2 = {:.4f}'.format(popt_c_power3[0] * popt_c_power3[0], R2_c_power3))
            #            plt.plot(x_new, biinversepower(x_new, *popt_c_bipower), '-', label='ec bi inverse power fitting: y = {:.4f} / x^{:.4f} + {:.4f} / x^{:.4f}  \n. R^2 = {:.4f}'.format(popt_c_bipower[0] * popt_c_bipower[0], popt_c_bipower[1], popt_c_bipower[2] * popt_c_bipower[2], popt_c_bipower[3] , R2_c_bipower))
            #             axs[1].xlim(0,10)
            axs[1].set_xlim(0,x_floor_plus1)
            axs[1].set_ylim(-7,1)
            #             plt.ylim(-0.1,1.2)
            axs[1].set_xlabel("$\chi$", fontsize = 15)
            axs[1].set_ylabel("log(Relative absolute error)", fontsize = 15)
            axs[1].legend(fontsize=12)


            # here we want to make two Y axis for the first plot to emphesize the convergence speed of exchange energy
            # The first curve is the data whose number of grid point over no +nv is less than 3, and the second curve is 
            # the data whose umber of grid point over no +nv is greater than 3.
            # you can set the x cutoff as you want
            x_cutoff = 2
            data_after_x_3 = data[data['ratio_grid_num'] >= x_cutoff]
            data_before_x_3 = data[data['ratio_grid_num'] < x_cutoff]
            x_before_3 = data_before_x_3['ratio_grid_num']
            x_after_3 = data_after_x_3['ratio_grid_num']

            y_x_before_3 = data_before_x_3['PerError_x']
            y_c_before_3 = data_before_x_3['PerError_c']

            y_x_after_3 = data_after_x_3['PerError_x']
            y_c_after_3 = data_after_x_3['PerError_c']

            ymax_after_3plus1 = (math.floor(max(y_x_after_3)*100) + 1)/100


            secaxs = axs[0].twinx()
            secaxs.tick_params(axis="y",direction="in", width = 1.5, length = 5, labelsize=15)
            axs[0].plot(x_before_3,  y_x_before_3  , label = "x", lw = 2.5)
            axs[0].plot(x_before_3, y_c_before_3  , label = "c", lw = 2.5)
            secaxs.plot(x_after_3,  y_x_after_3 , label = "x", lw = 2.5)
            secaxs.plot(x_after_3, y_c_after_3 , label = "c", lw = 2.5)
            axs[0].set_xlim(0,x_floor_plus1)
            axs[0].set_ylim(0,1)
            secaxs.set_ylim(0,ymax_after_3plus1)
            # axs[0].yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, decimals=None, symbol='%', is_latex=False))
            # secaxs.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, decimals=None, symbol='%', is_latex=False))
            axs[0].set_yticklabels(["{:,.0%}".format(y) for y in axs[0].get_yticks()], fontsize=15)
            secaxs.set_yticklabels(["{:,.0%}".format(y) for y in secaxs.get_yticks()], fontsize=15)
            axs[0].legend(fontsize=15)
            axs[0].set_xlabel("$\chi$", fontsize = 15)
            axs[0].set_ylabel("Percentage of relative absolute error", fontsize = 15)
            axs[0].axvline(x=x_cutoff, ymin=0, ymax=100, color= 'gray', linestyle = '--')
            figs.savefig( figname, dpi=600)


# save the fitting paramenter to a file
df_save.to_csv('./partial_energy_fitting/fitting_parameters.csv')
