#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
from scipy.spatial import distance
import math
import os
import sys
from bootstrap import iteration 
from histgram import hist

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.cm import ScalarMappable

#os.mkdir("./hello")
#os.chdir('graphs')
#dirr = ['mp2', 'mp3', 'ccsd']
#mett = ['dz', 'tz','adz']
#testt = ['svdthc', 'gridthc', 'thc']
#for i in dirr:
#    os.mkdir(i)
#    os.chdir(i)
#    for j in mett:
#        os.mkdir(j)
#        os.chdir(j)
#        for k in testt:
#            os.mkdir(k)
#        os.chdir('..')
#    os.chdir('..')
#os.chdir('..')
#os.chdir('..')


if (len(sys.argv) <= 3):
    print("Usage: python .py mol basis method ftype SizeofPairGridPoint")
    exit()
    


mol = sys.argv[1]
basis = sys.argv[2]
method = sys.argv[3]
Ftype = sys.argv[4]
size = int(sys.argv[5]) # the size of grid point you use, usually 
if basis == "dz/":
    basis = ""
#starting_points = int(sys.argv[3])
#nv_no_2 = int(sys.argv[4])
#percentage = int(round(starting_points/nv_no_2 * 100, 0))

if (basis == "adz/"):
    basis_str = "adz"
elif (basis == "rtz/"):
    basis_str = "rtz"
elif (basis == "tz/"):
    basis_str = "tz"
else:
    basis_str = "dz"
print("Chao 0", basis_str, method, Ftype)
# dira = dir
path = os.getcwd() 
filename = path + "/" + basis + "energyAnalysis.csv"
pvtFile = path +  "/" + basis + "grid/pvt.dat"
gridFile = path +  "/" + basis + "grid/grid.dat"
savecsvFile = path + "/bootstrapping_" + mol + "_" + basis_str + "_" + method + "_" + Ftype + ".csv" 
# filename = '/users/chaoy/work/THC/B3LYP-D3_def2_TZVP_Geom_6.0/%s/data/%s/energyAnalysis_all.csv' % dir % 
# filename = "/users/chaoy/work/THC/B3LYP-D3_def2_TZVP_Geom_6.0/{}/data/{}energyAnalysis.csv".format(mol, basis)
# pvtFile = "/users/chaoy/work/THC/B3LYP-D3_def2_TZVP_Geom_6.0/{}/data/{}grid/pvt.dat".format(mol, basis)
# gridFile = "/users/chaoy/work/THC/B3LYP-D3_def2_TZVP_Geom_6.0/{}/data/{}grid/grid.dat".format(mol, basis)
#savecsvFile = "/users/chaoy/work/THC/B3LYP-D3_def2_TZVP_Geom_6.0/{}/data/bootstrapping_{}_{}_{}_{}.csv".format(mol,mol,basis_str, method, Ftype)
print("Reading data from ", filename)
print("save datafrmae to csv file:", savecsvFile)
# filename = '~/work/alkanes/%d/results.csv' % dir

# read cvs file
df = pd.read_csv(filename, index_col = False)
df['Method'] = df['Method'].str.replace(" ", "") # remove the white space in the pandas columns
df['Ftype'] = df['Ftype'].str.replace(" ", "") # remove the white space in the pandas columns
#df['npairs'] = df['npairs'].str.replace(" ", "") # remove the white space in the pandas columns
#df['num_point_keep'] = df['num_point_keep'].str.replace(" ", "") # remove the white space in the pandas columns


#read in pvt
#pvtFile = open("data/grid/out/pvt.dat")
#pvtFile = open('./../%d/tests/%s/data/tz/grid/out/pvt.dat' % (mol,test))
pvtFile = open(pvtFile )
pvt = np.loadtxt(pvtFile, delimiter=" ")

                
#read in grid
#gridFile = open("grid/out/grid.dat")
gridFile = open(gridFile)
grid = np.loadtxt(gridFile, delimiter=" ")



Methods = df["Method"].unique()
Ftypes = df["Ftype"].unique()   
npairs = df["npairs"].unique()   
nums_point_keep = np.delete(df["num_point_keep"].unique(), 0)

print(df)
print(Methods)
print(Ftypes)
print(npairs)
print(nums_point_keep)
nv_no_2 = max(nums_point_keep)

# combine boot streaping and histgram together
# return the average of hist, the std of hist and the bins
def bsh(data, bins = 10, num_iteration = 5000,  weights = []):
    h_ave = []
    std_ave = []
    h, b = hist(data, bins = bins, weights = weights)

    for i in range(len(h)):
        temp = h[i] 
        ave, ave_min, ave_max, av_std = iteration(temp, num_iteration)
        h_ave.append(ave)
        std_ave.append(av_std)
    return np.array(h_ave), np.array(std_ave), np.array(b) 


# append a row to the end of a dataframe
def df_append(df1, ave, std, atom_dist_ave, bins, percentage ,energytype):
    row = [mol, method, basis_str, Ftype, npair, percentage, num_point_keep, energytype]
    stat_list = list(np.concatenate(( ave, std,  atom_dist_ave, bins), axis=None))
    row.extend(stat_list)
    df1 = df1.append(pd.DataFrame([row], columns=name), ignore_index=True)
    return df1

# set the name of new csv file
name = ['mol', 'method', 'basis', 'Ftype', 'npair', '# * (nv+no)', 'num_point_keep', 'energy_type',
       'ave1', 'ave2', 'ave3', 'ave4', 'ave5', 'ave6', 'ave7', 'ave8', 'ave9', 'ave10',
       'ave11', 'ave12', 'ave13', 'ave14', 'ave15', 'ave16', 'ave17', 'ave18', 'ave19', 'ave20',
       'ave21', 'ave22', 'ave23', 'ave24', 'ave25', 'ave26', 'ave27', 'ave28', 'ave29', 'ave30',
       'ave31', 'ave32', 'ave33', 'ave34', 'ave35', 'ave36', 'ave37', 'ave38', 'ave39', 'ave40',
       'ave41', 'ave42', 'ave43', 'ave44', 'ave45', 'ave46', 'ave47', 'ave48', 'ave49', 'ave50',
       'std1', 'std2', 'std3', 'std4', 'std5', 'std6', 'std7', 'std8', 'std9', 'std10',
       'std11', 'std12', 'std13', 'std14', 'std15', 'std16', 'std17', 'std18', 'std19', 'std20',
       'std21', 'std22', 'std23', 'std24', 'std25', 'std26', 'std27', 'std28', 'std29', 'std30',
       'std31', 'std32', 'std33', 'std34', 'std35', 'std36', 'std37', 'std38', 'std39', 'std40',
       'std41', 'std42', 'std43', 'std44', 'std45', 'std46', 'std47', 'std48', 'std49', 'std50',
       'ave_dist_to_atom_1', 'ave_dist_to_atom_2', 'ave_dist_to_atom_3', 'ave_dist_to_atom_4', 'ave_dist_to_atom_5','ave_dist_to_atom_6','ave_dist_to_atom_7','ave_dist_to_atom_8','ave_dist_to_atom_9','ave_dist_to_atom_10',
       'ave_dist_to_atom_11', 'ave_dist_to_atom_12', 'ave_dist_to_atom_13', 'ave_dist_to_atom_14', 'ave_dist_to_atom_15','ave_dist_to_atom_16','ave_dist_to_atom_17','ave_dist_to_atom_18','ave_dist_to_atom_19','ave_dist_to_atom_20',
       'ave_dist_to_atom_21', 'ave_dist_to_atom_22', 'ave_dist_to_atom_23', 'ave_dist_to_atom_24', 'ave_dist_to_atom_25','ave_dist_to_atom_26','ave_dist_to_atom_27','ave_dist_to_atom_28','ave_dist_to_atom_29','ave_dist_to_atom_30',
       'ave_dist_to_atom_31', 'ave_dist_to_atom_32', 'ave_dist_to_atom_33', 'ave_dist_to_atom_34', 'ave_dist_to_atom_35','ave_dist_to_atom_36','ave_dist_to_atom_37','ave_dist_to_atom_38','ave_dist_to_atom_39','ave_dist_to_atom_40',
       'ave_dist_to_atom_41', 'ave_dist_to_atom_42', 'ave_dist_to_atom_43', 'ave_dist_to_atom_44', 'ave_dist_to_atom_45','ave_dist_to_atom_46','ave_dist_to_atom_47','ave_dist_to_atom_48','ave_dist_to_atom_49','ave_dist_to_atom_50',
       'bin1', 'bin2', 'bin3', 'bin4', 'bin5', 'bin6', 'bin7', 'bin8', 'bin9', 'bin10',
       'bin11', 'bin12', 'bin13', 'bin14', 'bin15', 'bin16', 'bin17', 'bin18', 'bin19', 'bin20',
       'bin21', 'bin22', 'bin23', 'bin24', 'bin25', 'bin26', 'bin27', 'bin28', 'bin29', 'bin30',
       'bin31', 'bin32', 'bin33', 'bin34', 'bin35', 'bin36', 'bin37', 'bin38', 'bin39', 'bin40',
       'bin41', 'bin42', 'bin43', 'bin44', 'bin45', 'bin46', 'bin47', 'bin48', 'bin49', 'bin50', 'bin51']

df_new = pd.DataFrame(columns=name)


bootstreapinginterations = 5000
print(nums_point_keep)
print(nv_no_2)
for npair in npairs:
    for num_point_keep in nums_point_keep:
        #data = df[(df['Method'] == method) & df['Ftype'] == Ftype & df['npairs'] == npair & (df['num_point_keep'] == num_point_keep)]
        data = df[(df['Method'] == method) & (df['Ftype'] == Ftype) &  (df['npairs'] == npair) & (df['num_point_keep'] == num_point_keep)] 
        percentage = round(num_point_keep/nv_no_2 * 6, 1)
        print(method, Ftype, npair, num_point_keep, percentage)
        # data2 = data
        
        #if not os.path.exists('./../%d/tests/eS_mp2_tz/graph' % (dir)):
            #os.mkdir('./../%d/tests/eS_mp2_tz/graph' % (dir))
        
        
        #print(datav)
        #dz - 582
        #adz - 1014
        #tz - 1452
        
        #print(data) #make a plot as increase number of single points, plot of thc c vs exact c and x vs exact x
        
        #create the numpy matrices for v and c
        v4x = data[['Point1', 'Point2', 'E4X']]
        v4c = data[['Point1', 'Point2', 'E4C']]
        v8x = data[['Point1', 'Point2', 'E8X']]
        v8c = data[['Point1', 'Point2', 'E8C']]
        fval = data['Fval']
        
        v2c = data['E2C']
        v2x = data['E2X']
        
        v4xdata = v4x.to_numpy()
        v4cdata = v4c.to_numpy()
        v8xdata = v8x.to_numpy()
        v8cdata = v8c.to_numpy()
        fvaldata = fval.to_numpy()
        
        v4xdata[:,2] = v4xdata[:,2]/v2x
        v4cdata[:,2] = v4cdata[:,2]/v2c
        v8xdata[:,2] = v8xdata[:,2]/v2x
        v8cdata[:,2] = v8cdata[:,2]/v2c
        
        v4data = v4x.to_numpy()
        v8data = v8x.to_numpy()
        vxdata = v4x.to_numpy()
        vcdata = v4c.to_numpy()
        vx_vdata = v4x.to_numpy()
        vc_vdata = v4c.to_numpy()
        vdata = v4x.to_numpy()
        vx_vcdata = v4x.to_numpy()
        
        v4data[:,2] = (v4xdata[:,2]*v2x + v4cdata[:,2]*v2c)/(v2c + v2x)
        v8data[:,2] = (v8xdata[:,2]*v2x + v8cdata[:,2]*v2c)/(v2c + v2x)
        vxdata[:,2] = (v8xdata[:,2] + v4xdata[:,2])
        vcdata[:,2] = (v8cdata[:,2] + v4cdata[:,2])
        vx_vdata[:,2] = (v4xdata[:,2]*v2x + v8xdata[:,2]*v2x) / (v2c + v2x)
        vc_vdata[:,2] = (v4cdata[:,2]*v2c + v8cdata[:,2]*v2c) / (v2c + v2x)
        vdata[:,2] = (v4xdata[:,2]*v2x + v4cdata[:,2]*v2c + v8xdata[:,2]*v2x + v8cdata[:,2]*v2c)/(v2c + v2x)
        vx_vcdata[:,2] = vxdata[:,2] - vcdata[:,2]
        

        # v4x = data.drop(columns = ['Exact_c','Exact_x', 'LS-THC_c','LS-THC_x','E4_c','E8_c', 'E8_x','S'])
        # v4c = data2.drop(columns = ['Exact_c','Exact_x', 'LS-THC_c','LS-THC_x','E4_x','E8_c', 'E8_x','S'])
        # v8c = data.drop(columns = ['Exact_c','Exact_x', 'LS-THC_c','LS-THC_x','E4_c','E4_x', 'E8_x','S'])
        # v8x = data2.drop(columns = ['Exact_c','Exact_x', 'LS-THC_c','LS-THC_x','E4_x','E4_c','E8_c','S'])
        # #c = datac.drop(columns = ['target', 'basis', 'calc', 'test', 'rank', 'P'])
        # #v_single = datav_single.drop(columns = ['target', 'basis', 'calc', 'test', 'rank', 'P'])
        # #c_single = datac_single.drop(columns = ['target', 'basis', 'calc', 'test', 'rank','P'])
        
        # v4xdata = v4x.to_numpy()
        # v4cdata = v4c.to_numpy()
        # v8xdata = v8x.to_numpy()
        # v8cdata = v8c.to_numpy()
        
        # v4data = v4x.to_numpy()
        # v8data = v8x.to_numpy()
        # vxdata = v4x.to_numpy()
        # vcdata = v4c.to_numpy()
        # vdata = v4x.to_numpy()
        
        # v4data[:,2] = v4xdata[:,2] + v4cdata[:,2]
        # v8data[:,2] = v8xdata[:,2] + v8cdata[:,2]
        # vxdata[:,2] = v8xdata[:,2] + v4xdata[:,2]
        # vcdata[:,2] = v8cdata[:,2] + v4cdata[:,2]
        # vdata[:,2] = vcdata[:,2] + vxdata[:,2]
        
        # test = 'energyS_ccsd_tz'
        
        # #cdata = c.to_numpy()
        # #vdata_single = v_single.to_numpy()
        # #cdata_single = c_single.to_numpy()
        # size = np.shape(v4data)[0]
        # #size_s = np.shape(vdata_single)[0]
        
        
        
        
        #make the lookup matrices
        lookup = np.zeros((size, 2)) #2000 tests, 2 points each
        
        for i in range(size):
            lookup[i][0] = pvt[int(v4xdata[i][0])]
            lookup[i][1] = pvt[int(v4xdata[i][1])]
        
        
        
        
        print("Chao 4")
        #obtain x,y,z and weight for each point in each pair
        point1 = np.zeros((size,3))
        point2 = np.zeros((size,3))
        weight1 = np.zeros(size)
        weight2 = np.zeros(size)
        
        
        for i in range(size):
            point1[i][:] = grid[int(lookup[i][0])][0:3]
            point2[i][:] = grid[int(lookup[i][1])][0:3]
            weight1[i] = grid[int(lookup[i][0])][3]
            weight2[i] = grid[int(lookup[i][1])][3]
        
        print("Chao 3")
        print(point1)
        print(point2)
        
        
        
        #distance matrix
        dist = distance.cdist(point1, point2, 'euclidean')
        # dir = 8
        
        print("Chao 5")
        
        
        def GetAtomLoc(ZMATfilename):
            loclist = []
            with open(ZMATfilename) as file:
                next(file)
                for line in file:
                    if not line.strip(): #if line is blank. 
                        break
                    loclist.append(line.split())
            # get rid of the first column of list, which is the atom name
            loclisttmp = np.float_([(row[1:]) for row in loclist])
            return loclisttmp
        
        
        # dir = dira
        # dir = "adenine"
        atomFile = path + "/../ZMAT"
        # atomFile = "/users/chaoy/work/THC/B3LYP-D3_def2_TZVP_Geom_6.0/%s/ZMAT"%(mol)
        atom = GetAtomLoc(atomFile)
        # atom = np.loadtxt(atomFile, delimiter=" ")
        
        # create distance matrices from atoms
        print("Calculting distance")
        atomDist1 = distance.cdist(point1, atom, 'euclidean')
        atomDist2 = distance.cdist(point2, atom, 'euclidean')
        print("Calculting Done")
        
        # get min along atom axis and average
        min1 = np.min(atomDist1, axis=1)
        min2 = np.min(atomDist2, axis=1)
        
        avg = (min1 + min2) / 2
        # dir = 8
        
        
        
        # make the folder for Dis Vs Energy
        # './pic/%s/%s/%s/DisVSEx/ExvsEMP2x_%s_%s_%s_%s_%i.png' %(basis_str, method, Ftype, mol,  basis_str,     method, Ftype, percentage))
        dis_vs_ene_path = "./pic/%s/%s/%s/EaddoverE"%(basis_str, method, Ftype)
        if (not os.path.exists(dis_vs_ene_path)):
            os.makedirs(dis_vs_ene_path)
        
        dis_vs_fval_path = "./pic/%s/%s/%s/Fval"%(basis_str, method, Ftype)
        if (not os.path.exists(dis_vs_fval_path)):
            os.makedirs(dis_vs_fval_path)
        
        dis_vs_ex_path = "./pic/%s/%s/%s/ExaddoverEx"%(basis_str, method, Ftype)
        if (not os.path.exists(dis_vs_ex_path)):
            os.makedirs(dis_vs_ex_path)
        
        dis_vs_ec_path = "./pic/%s/%s/%s/EcaddoverEc"%(basis_str, method, Ftype)
        if (not os.path.exists(dis_vs_ec_path)):
            os.makedirs(dis_vs_ec_path)

        dis_vs_exaddovere_path = "./pic/%s/%s/%s/ExaddoverE"%(basis_str, method, Ftype)
        if (not os.path.exists(dis_vs_exaddovere_path)):
            os.makedirs(dis_vs_exaddovere_path)

        dis_vs_ecaddovere_path = "./pic/%s/%s/%s/EcaddoverE"%(basis_str, method, Ftype)
        if (not os.path.exists(dis_vs_ecaddovere_path)):
            os.makedirs(dis_vs_ecaddovere_path)
        print("Chao 2")
        # weighted histogram of distances for e4x
        pairDistance = np.diag(dist)
        
        # weighted histogram of distances for e4x
        pairDistance = np.diag(dist)
        
        # ratio of distances for ex_add/ex_0
        # ratio of distances for ex
        h, std, b = bsh(pairDistance, bins = 50, num_iteration = bootstreapinginterations,  weights = vxdata[:, 2])
        hbar, b = hist(pairDistance, bins = 50, weights = avg, ave= True)
        h_m_std = h - std
        h_p_std = h + std
        scilim = math.floor(np.log10(max(abs(h)))) # get the limitation of sci notation of y axis
        hbarmax = hbar.max()
        # get the max or min value of hbar with the std
        ylim = max(h_p_std) * 1.1
        
        # update the dataframe
        df_new = df_append(df_new, h, std, hbar, b, percentage, 'Ex_add/Ex')
        
        cm = plt.cm.get_cmap('RdYlBu_r')
        # cm = plt.cm.get_cmap('Set1')
        
        # Get the histogramp
        hbar_span = math.floor(hbarmax) + 1
        C = ([cm(((x)/hbar_span)) for x in hbar])
        sm = plt.cm.ScalarMappable(cmap=cm)
        ticks = np.linspace(0, hbar_span, 6, endpoint=True)
        
        fig, ax = plt.subplots(figsize=(7, 3.5), dpi = 600)
        plt.rcParams['axes.linewidth'] = 1
        _, _, bars = ax.hist(b[:-1], bins= b, weights = h, rwidth = 0.6)
        ax.fill_between(b[:-1]+(b[1] - b[0])/2, h_p_std, h_m_std, alpha=0.2, color = 'gray')
        
        # set the corlorbar scheme
        for bar, color in zip(bars, C):
            bar.set_facecolor(color)
        cbar = fig.colorbar(sm, orientation='vertical', fraction = 0.03)
        cbar.set_ticklabels(np.round((ticks), 2), fontsize = 10)
        
        # set the x y axis labels and ticks
        ax.set_xlabel("Distance between pairs of points (a.u.)", fontsize = 10, fontfamily = 'sans-serif')
        ax.ticklabel_format(axis='y', style='sci', scilimits = (scilim,scilim))
        ax.tick_params(axis="y",direction="in", width = 1, length = 3, labelsize=10)
        ax.tick_params(axis="x",direction="in", width = 1, length = 3, labelsize=10)

        ax.set_ylim(0, ylim)
        fig.text(0.45,0.8, "$\chi$= {}".format(percentage), fontfamily = 'sans-serif', fontweight ='extra bold', )
        ax.set_ylabel("$\\frac{\Delta E_x}{E_x}$", rotation = 'horizontal', fontsize = 10, fontfamily = 'sans-serif', labelpad = 5, )

        ax.set_xlim(0,max(b[:-1]) * 1.1) # set the x lim
        plt.savefig('./pic/%s/%s/%s/ExaddoverEx/ExaddoverEx_%s_%s_%s_%s_%s_%.1f.png' %(basis_str, method, Ftype, method, mol, basis_str,method, Ftype, percentage))
        
        
        # ratio of distances for ec_add/ec_0
        h, std, b = bsh(pairDistance, bins = 50, num_iteration = bootstreapinginterations,  weights = vcdata[:, 2])
        hbar, b = hist(pairDistance, bins = 50, weights = avg, ave= True)
        h_m_std = h - std
        h_p_std = h + std
        scilim = math.floor(np.log10(max(abs(h)))) # get the limitation of sci notation of y axis
        hbarmax = hbar.max()
        # get the max or min value of hbar with the std
        ylim = max(h_p_std) * 1.1
        
        # update the dataframe
        df_new = df_append(df_new, h, std, hbar, b,  percentage, 'Ec_add/Ec')
        
        cm = plt.cm.get_cmap('RdYlBu_r')
        # cm = plt.cm.get_cmap('Set1')
        
        # Get the histogramp
        hbar_span = math.floor(hbarmax) + 1
        C = ([cm(((x)/hbar_span)) for x in hbar])
        sm = plt.cm.ScalarMappable(cmap=cm)
        ticks = np.linspace(0, hbar_span, 6, endpoint=True)
        
        fig, ax = plt.subplots(figsize=(7, 3.5), dpi = 600)
        plt.rcParams['axes.linewidth'] = 1
        _, _, bars = ax.hist(b[:-1], bins= b, weights = h, rwidth = 0.6)
        ax.fill_between(b[:-1]+(b[1] - b[0])/2, h_p_std, h_m_std, alpha=0.2, color = 'gray')
        
        # set the corlorbar scheme
        for bar, color in zip(bars, C):
            bar.set_facecolor(color)
        cbar = fig.colorbar(sm, orientation='vertical', fraction = 0.03)
        cbar.set_ticklabels(np.round((ticks), 2), fontsize = 10)
        
        # set the x y axis labels and ticks
        ax.set_xlabel("Distance between pairs of points (a.u.)", fontsize = 10, fontfamily = 'sans-serif')
        ax.ticklabel_format(axis='y', style='sci', scilimits = (scilim,scilim))
        ax.tick_params(axis="y",direction="in", width = 1, length = 3, labelsize=10)
        ax.tick_params(axis="x",direction="in", width = 1, length = 3, labelsize=10)

        ax.set_ylim(0, ylim)
        fig.text(0.45,0.8, "$\chi$= {}".format(percentage), fontfamily = 'sans-serif', fontweight ='extra bold', )
        ax.set_ylabel("$\\frac{\Delta E_c}{E_c}$", rotation = 'horizontal', fontsize = 10, fontfamily = 'sans-serif', labelpad = 5, )

        ax.set_xlim(0,max(b[:-1]) * 1.1) # set the x lim
        plt.savefig('./pic/%s/%s/%s/EcaddoverEc/EcaddvsEc_%s_%s_%s_%s_%s_%.1f.png' %(basis_str, method, Ftype, method, mol, basis_str,method, Ftype, percentage))
        
        
        # ratio of distances for ex_add/e_0  (ex_0+ec_0)
        h, std, b = bsh(pairDistance, bins = 50, num_iteration = bootstreapinginterations,  weights = vx_vdata[:, 2])
        hbar, b = hist(pairDistance, bins = 50, weights = avg, ave= True)
        h_m_std = h - std
        h_p_std = h + std
        scilim = math.floor(np.log10(max(abs(h)))) # get the limitation of sci notation of y axis
        hbarmax = hbar.max()
        # get the max or min value of hbar with the std
        ylim = min(h_m_std) * 1.1
        
        # update the dataframe
        df_new = df_append(df_new, h, std, hbar, b, percentage, 'Ex_add/E')

        cm = plt.cm.get_cmap('RdYlBu_r')
        # cm = plt.cm.get_cmap('Set1')
        
        # Get the histogramp
        hbar_span = math.floor(hbarmax) + 1
        C = ([cm(((x)/hbar_span)) for x in hbar])
        sm = plt.cm.ScalarMappable(cmap=cm)
        ticks = np.linspace(0, hbar_span, 6, endpoint=True)
        
        fig, ax = plt.subplots(figsize=(7, 3.5), dpi = 600)
        plt.rcParams['axes.linewidth'] = 1
        _, _, bars = ax.hist(b[:-1], bins= b, weights = h, rwidth = 0.6)
        ax.fill_between(b[:-1]+(b[1] - b[0])/2, h_p_std, h_m_std, alpha=0.2, color = 'gray')
        
        # set the corlorbar scheme
        for bar, color in zip(bars, C):
            bar.set_facecolor(color)
        cbar = fig.colorbar(sm, orientation='vertical', fraction = 0.03)
        cbar.set_ticklabels(np.round((ticks), 2), fontsize = 10)
        
        # set the x y axis labels and ticks
        ax.set_xlabel("Distance between pairs of points (a.u.)", fontsize = 10, fontfamily = 'sans-serif')
        ax.ticklabel_format(axis='y', style='sci', scilimits = (scilim,scilim))
        ax.tick_params(axis="y",direction="in", width = 1, length = 3, labelsize=10)
        ax.tick_params(axis="x",direction="in", width = 1, length = 3, labelsize=10)

        ax.set_ylim(ylim, 0)
        fig.text(0.45,0.8, "$\chi$= {}".format(percentage), fontfamily = 'sans-serif', fontweight ='extra bold', )
        ax.set_ylabel("$\\frac{\Delta E_x}{E}$", rotation = 'horizontal', fontsize = 10, fontfamily = 'sans-serif', labelpad = 5, )

        ax.set_xlim(0,max(b[:-1]) * 1.1) # set the x lim
        
        plt.savefig('./pic/%s/%s/%s/ExaddoverE/ExaddvsE_%s_%s_%s_%s_%s_%.1f.png' %(basis_str, method, Ftype, method, mol, basis_str,method, Ftype, percentage))
        
        
        
        # ratio of distances for ec_add/(ec_0 + ex_0)
        h, std, b = bsh(pairDistance, bins = 50, num_iteration = bootstreapinginterations,  weights = vc_vdata[:, 2])
        hbar, b = hist(pairDistance, bins = 50, weights = avg, ave= True)
        h_m_std = h - std
        h_p_std = h + std
        scilim = math.floor(np.log10(max(abs(h)))) # get the limitation of sci notation of y axis
        hbarmax = hbar.max()
        # get the max or min value of hbar with the std
        ylim = max(h_p_std) * 1.1
        
        
        df_new = df_append(df_new, h, std, hbar, b,  percentage, 'Ec_add/E')
        cm = plt.cm.get_cmap('RdYlBu_r')
        # cm = plt.cm.get_cmap('Set1')
        
        # Get the histogramp
        hbar_span = math.floor(hbarmax) + 1
        C = ([cm(((x)/hbar_span)) for x in hbar])
        sm = plt.cm.ScalarMappable(cmap=cm)
        ticks = np.linspace(0, hbar_span, 6, endpoint=True)
        
        fig, ax = plt.subplots(figsize=(7, 3.5), dpi = 600)
        plt.rcParams['axes.linewidth'] = 1
        _, _, bars = ax.hist(b[:-1], bins= b, weights = h, rwidth = 0.6)
        ax.fill_between(b[:-1]+(b[1] - b[0])/2, h_p_std, h_m_std, alpha=0.2, color = 'gray')
        
        # set the corlorbar scheme
        for bar, color in zip(bars, C):
            bar.set_facecolor(color)
        cbar = fig.colorbar(sm, orientation='vertical', fraction = 0.03)
        cbar.set_ticklabels(np.round((ticks), 2), fontsize = 10)
        
        # set the x y axis labels and ticks
        ax.set_xlabel("Distance between pairs of points (a.u.)", fontsize = 10, fontfamily = 'sans-serif')
        ax.ticklabel_format(axis='y', style='sci', scilimits = (scilim,scilim))
        ax.tick_params(axis="y",direction="in", width = 1, length = 3, labelsize=10)
        ax.tick_params(axis="x",direction="in", width = 1, length = 3, labelsize=10)

        ax.set_ylim(0, ylim)
        fig.text(0.45,0.8, "$\chi$= {}".format(percentage), fontfamily = 'sans-serif', fontweight ='extra bold', )
        ax.set_ylabel("$\\frac{\Delta E_c}{E}$", rotation = 'horizontal', fontsize = 10, fontfamily = 'sans-serif', labelpad = 5, )

        ax.set_xlim(0,max(b[:-1]) * 1.1) # set the x lim
        
        plt.savefig('./pic/%s/%s/%s/EcaddoverE/EcaddoverE_%s_%s_%s_%s_%s_%.1f.png' %(basis_str, method, Ftype, method, mol, basis_str,method, Ftype, percentage))
        
        # ratio of distances for total
        h, std, b = bsh(pairDistance, bins = 50, num_iteration = bootstreapinginterations,  weights = vdata[:, 2])
        hbar, b = hist(pairDistance, bins = 50, weights = avg, ave= True)
        h_m_std = h - std
        h_p_std = h + std
        scilim = math.floor(np.log10(max(abs(h)))) # get the limitation of sci notation of y axis
        hbarmax = hbar.max()
        # get the max or min value of hbar with the std
        ymaxlim = max(h_p_std) * 1.1
        yminlim = min(h_m_std) * 1.1
        
        df_new = df_append(df_new, h, std, hbar, b,  percentage, 'E_add/E')
        cm = plt.cm.get_cmap('RdYlBu_r')
        # cm = plt.cm.get_cmap('Set1')
        
        # Get the histogramp
        hbar_span = math.floor(hbarmax) + 1
        C = ([cm(((x)/hbar_span)) for x in hbar])
        sm = plt.cm.ScalarMappable(cmap=cm)
        ticks = np.linspace(0, hbar_span, 6, endpoint=True)
        
        fig, ax = plt.subplots(figsize=(7, 3.5), dpi = 600)
        plt.rcParams['axes.linewidth'] = 1
        _, _, bars = ax.hist(b[:-1], bins= b, weights = h, rwidth = 0.6)
        ax.fill_between(b[:-1]+(b[1] - b[0])/2, h_p_std, h_m_std, alpha=0.2, color = 'gray')
        
        # set the corlorbar scheme
        for bar, color in zip(bars, C):
            bar.set_facecolor(color)
        cbar = fig.colorbar(sm, orientation='vertical', fraction = 0.03)
        cbar.set_ticklabels(np.round((ticks), 2), fontsize = 10)
        
        # set the x y axis labels and ticks
        ax.set_xlabel("Distance between pairs of points (a.u.)", fontsize = 10, fontfamily = 'sans-serif')
        ax.ticklabel_format(axis='y', style='sci', scilimits = (scilim,scilim))
        ax.tick_params(axis="y",direction="in", width = 1, length = 3, labelsize=10)
        ax.tick_params(axis="x",direction="in", width = 1, length = 3, labelsize=10)
        if (ymaxlim > 0 and yminlim > 0):
            ax.set_ylim(0, ymaxlim)
        elif (ymaxlim < 0 and yminlim < 0):
            ax.set_ylim(yminlim, 0)
        else:
            ax.set_ylim(yminlim, ymaxlim)
        fig.text(0.45,0.8, "$\chi$= {}".format(percentage), fontfamily = 'sans-serif', fontweight ='extra bold', )
        ax.set_ylabel("$\\frac{\Delta E}{E}$", rotation = 'horizontal', fontsize = 10, fontfamily = 'sans-serif', labelpad = 5, )

        ax.set_xlim(0,max(b[:-1]) * 1.1) # set the x lim
        
        plt.savefig('./pic/%s/%s/%s/EaddoverE/EaddoverE_%s_%s_%s_%s_%s_%.1f.png' %(basis_str, method, Ftype, method, mol,  basis_str, method, Ftype, percentage))
        
        
        
        # ratio of distances for Fval
        h, std, b = bsh(pairDistance, bins = 50, num_iteration = bootstreapinginterations,  weights = fvaldata)
        hbar, b = hist(pairDistance, bins = 50, weights = avg, ave= True)
        h_m_std = h - std
        h_p_std = h + std
        scilim = math.floor(np.log10(max(abs(h)))) # get the limitation of sci notation of y axis
        hbarmax = hbar.max()
        # get the max or min value of hbar with the std
        ylim = max(h_p_std) * 1.1
        
        df_new = df_append(df_new, h, std, hbar, b,  percentage, 'Fval')
        cm = plt.cm.get_cmap('RdYlBu_r')
        # cm = plt.cm.get_cmap('Set1')
        
        # Get the histogramp
        hbar_span = math.floor(hbarmax) + 1
        C = ([cm(((x)/hbar_span)) for x in hbar])
        sm = plt.cm.ScalarMappable(cmap=cm)
        ticks = np.linspace(0, hbar_span, 6, endpoint=True)
        
        fig, ax = plt.subplots(figsize=(7, 3.5), dpi = 600)
        plt.rcParams['axes.linewidth'] = 1
        _, _, bars = ax.hist(b[:-1], bins= b, weights = h, rwidth = 0.6)
        ax.fill_between(b[:-1]+(b[1] - b[0])/2, h_p_std, h_m_std, alpha=0.2, color = 'gray')
        
        # set the corlorbar scheme
        for bar, color in zip(bars, C):
            bar.set_facecolor(color)
        cbar = fig.colorbar(sm, orientation='vertical', fraction = 0.03)
        cbar.set_ticklabels(np.round((ticks), 2), fontsize = 10)
        
        # set the x y axis labels and ticks
        ax.set_xlabel("Distance between pairs of points (a.u.)", fontsize = 10, fontfamily = 'sans-serif')
        ax.ticklabel_format(axis='y', style='sci', scilimits = (scilim,scilim))
        ax.tick_params(axis="y",direction="in", width = 1, length = 3, labelsize=10)
        ax.tick_params(axis="x",direction="in", width = 1, length = 3, labelsize=10)

        ax.set_ylim(0, ylim)
        fig.text(0.45,0.8, "$\chi$= {}".format(percentage), fontfamily = 'sans-serif', fontweight ='extra bold', )
        ax.set_ylabel("$f$", rotation = 'horizontal', fontsize = 15, fontfamily = 'sans-serif', labelpad = 5, )

        ax.set_xlim(0,max(b[:-1]) * 1.1) # set the x lim
        
        plt.savefig('./pic/%s/%s/%s/Fval/Fval_%s_%s_%s_%s_%.1f.png' %(basis_str, method, Ftype, mol,  basis_str,method, Ftype, percentage))
        
        #plt.savefig(f'graphs/{method}/{basis}/{test}/ratio_distance_c_{p}.png')
        # plt.savefig('./../%d/tests/%s/graph/ratio_distance_total.png' %(mol,test))
        
        
        
        # unweighted geometric mean of weights
        geoWeights = np.log10(np.sqrt(abs(np.multiply(weight1, weight2))))
        
        
        # get min along atom axis and average
        min1 = np.min(atomDist1, axis=1)
        min2 = np.min(atomDist2, axis=1)
        
        avg = (min1 + min2) / 2
        
        
        # make the folder for hist of the ave distance of pairs of points to its cloest atom 
        ave_vs_ene_path = "./pic/%s/%s/%s/Avedist"%(basis_str, method, Ftype)
        if (not os.path.exists(ave_vs_ene_path)):
            os.makedirs(ave_vs_ene_path)
        
        
        # histogram of distances from atoms for v  and c
        fig, ax = plt.subplots(figsize=(10, 6), dpi =600)
        # ax.set_title('Unweighted Histogram of Average Distances from Atoms ({})'.format(basis_str), fontsize = 17)
        ax.hist(avg, rwidth = 0.5, bins = 50)
        ax.set_ylabel("#count", fontsize = 15)
        ax.set_xlabel("Average distances between pair of grid points to its closest atom (a.)", fontsize = 10)
        #plt.title('Unweighted Histogram of Average Distances from Atoms ({})'.format(basis_str), fontsize = 20)
        #plt.savefig(f'graphs/{method}/{basis}/{test}/unweighted_avgdist_{p}.png')
        plt.savefig('./pic/%s/%s/%s/Avedist/Avedist_%s_%s_%s_%s.png' %(basis_str, method, Ftype, mol,  basis_str, method, Ftype))
        #plt.savefig('./pic/%s/%s/%s/DisVSEx/ExvsEMP2x_%s_%s_%s_%s_%i.png' %(basis_str, method, Ftype, mol,  basis_str,method, Ftype, percentage))
        
        
        
df_new.to_csv(savecsvFile)
