'''
Bootstraping in python. created by Chao Yin in Mar 20, 2023
Prerequisite: numpy (you can use canda or pip3 to install it)
Usage: put the dataset you have and the number of iterations you want
The return is some statistics values (ave, min, max, std) of the bootstraping.
'''

from random import uniform
from numpy import std

def sampling(init_sample, num_s = 1000):
    sample = []
    size_s = len(init_sample)
    for i in range(size_s):
        d = init_sample[int(uniform(0,size_s))]
        sample.append(d)
    return sum(sample)/size_s

def iteration(init_sample, num_s = 1000):
    size_s = len(init_sample)
    if size_s == 0:
        return 0.0, 0.0, 0.0, 0.0
    else:
        ave_list = []
        for j in range(num_s):
            ave = sampling(init_sample, num_s)
            ave_list.append(ave)
        ave_smapling = sum(ave_list)/num_s
        ave_min = min(ave_list)
        ave_max = max(ave_list)
        ave_std = std(ave_list)
        return ave_smapling, ave_min, ave_max, ave_std


if __name__ == "__main__":
    data = [1,2,2,3,3,3,4,4,4,4]
    num_it = 1000
    ave, ave_min, ave_max, av_std = iteration(data, num_it)
    print("ave, ave_min, ave_max, av_std = ", ave, ave_min, ave_max, av_std)
