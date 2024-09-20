import numpy as np


def hist(a, bins = 10, weights = [], ave = False):
    amax = max(a)
    amin = min(a)
    intverals = np.linspace(amin, amax+0.001, num = bins+ 1)
    inds = np.digitize(a, intverals)

    if (len(weights) == 0):
        unweighted_list =  [[] for i in range(bins)]
        for i in range(len(a)):
            ind = inds[i]
            unweighted_list[ind-1].append(1)
        # if the ave option is False, it will return the invidival values for the histgram, otherwise it will return the average vale.
        if (ave == False):
            return unweighted_list,  intverals
        else:
            unweighted_list_ave = []
            for j in range(bins):
                if (len(unweighted_list[j]) == 0):
                    unweighted_list_ave.append(0.0)
                else:
                    ave_temp = np.mean(unweighted_list[j])
                    unweighted_list_ave.append(ave_temp)
            return np.array(unweighted_list_ave), intverals

    else:
        weighted_list =  [[] for i in range(bins)]
        for i in range(len(a)):
            ind = inds[i]
            weighted_list[ind-1].append(1 * weights[i])
        if (ave == False):
            return weighted_list,  intverals
        else:
            weighted_list_ave = []
            for j in range(bins):
                if (len(weighted_list[j]) == 0):
                    weighted_list_ave.append(0.0)
                else:
                    ave_temp = np.mean(weighted_list[j])
                    weighted_list_ave.append(ave_temp)
            return np.array(weighted_list_ave), intverals



if __name__ == "__main__":
    data = np.array([1,2,3,4,5,6,7,8,9,10])
    weights = np.array([0.1, 0.1, 0.3, 0.2, 0.1, 0.5,0.1, 0.2, 0.1, 0.2])
    res, intevrals = hist(data, bins = 2) #, weights = weights)
    print(res)
    print(intevrals)


    
            

                


