from struct import unpack
import os
from shutil import copy
import numpy as np


os.mkdir("data")


types = ["adz",".","tz","rtz"]
dirs = [["AUG-PVDZ-CCSD","AUG-PVDZ-MP2","AUG-PVDZ-MP3"],["PVDZ-CCSD","PVDZ-MP2","PVDZ-MP3"],["PVTZ-CCSD","PVTZ-MP2","PVTZ-MP3"],["REDUCED_PVTZ-CCSD","REDUCED_PVTZ-MP2","REDUCED_PVTZ-MP3"]]


for i in range(4):
    os.chdir("data")
    type = types[i]
    exit = "./.."
    if type != '.':
        os.mkdir(type)
        os.chdir(type)
        exit = "./../.."
    
    os.mkdir("grid")
    os.mkdir("ccsd")

    copy(('%s/%s/fa.dat')%(exit, dirs[i][2]),".")
    copy(('%s/%s/fi.dat')%(exit, dirs[i][2]),".")
    copy(('%s/%s/v.dat')%(exit, dirs[i][2]),".")
    copy(('%s/%s/t.dat')%(exit, dirs[i][2]),".")
    os.chdir("grid")
    copy(('%s/../%s/XA.dat')%(exit, dirs[i][1]),".")
    copy(('%s/../%s/XI.dat')%(exit, dirs[i][1]),".")
    copy(('%s/../%s/grid.dat')%(exit, dirs[i][1]),".")
    copy(('%s/../%s/pvt.dat')%(exit, dirs[i][1]),".")
    os.chdir("../ccsd")
    copy(('%s/../%s/fa.dat')%(exit, dirs[i][0]),".")
    copy(('%s/../%s/fi.dat')%(exit, dirs[i][0]),".")
    copy(('%s/../%s/v.dat')%(exit, dirs[i][0]),".")
    copy(('%s/../%s/t.dat')%(exit, dirs[i][0]),".")
    os.chdir("..")

    os.chdir(exit)
'''   
for i in range(4):
    os.chdir("data")
    type = types[i]
    exit = "./.."
    if type != '.':
        os.chdir(type)
        exit = "./../.."

    os.remove(('%s/%s/fa.dat')%(exit, dirs[i][2]))
    os.remove(('%s/%s/fi.dat')%(exit, dirs[i][2]))
    os.remove(('%s/%s/v.dat')%(exit, dirs[i][2]))
    os.remove(('%s/%s/t.dat')%(exit, dirs[i][2]))
    os.chdir("grid")
    os.remove(('%s/../%s/XA.dat')%(exit, dirs[i][1]))
    os.remove(('%s/../%s/XI.dat')%(exit, dirs[i][1]))
    os.remove(('%s/../%s/grid.dat')%(exit, dirs[i][1]))
    os.remove(('%s/../%s/pvt.dat')%(exit, dirs[i][1]))
    os.chdir("../ccsd")
    os.remove(('%s/../%s/fa.dat')%(exit, dirs[i][0]))
    os.remove(('%s/../%s/fi.dat')%(exit, dirs[i][0]))
    os.remove(('%s/../%s/v.dat')%(exit, dirs[i][0]))
    os.remove(('%s/../%s/t.dat')%(exit, dirs[i][0]))
    os.chdir("..")

    os.chdir(exit)
'''  