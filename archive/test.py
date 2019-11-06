import os, numpy, csv, operator
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from numpy import mean,std
from sklearn.feature_selection import VarianceThreshold
from sklearn.preprocessing import StandardScaler
from sklearn.externals import joblib as pickle
from sklearn.decomposition import PCA,SparsePCA,MiniBatchSparsePCA,TruncatedSVD
from sklearn.decomposition import RandomizedPCA as RCA
# from paramRegressor import *
from collections import OrderedDict

def loadData(name,path='.'):
    '''
    This loads a pickle file and returns the content which is a DICTIONARY object in our case.
    '''
    if ".pkl" in name:
            name = name.split(".pkl")[0]
    if "/" in name:
            name = name.split("/",1)[1]

    return pickle.load(path+"/"+name + '.pkl')
    # with open(path+"/"+name + '.pkl', 'rb') as f:
    #       return pickle.load(f)

def saveData(obj, name,path='.'):
    '''
    This saves a object into a pickle file. In our case, it is generally a DICTIONARY object.
    '''

    with open(path+"/"+name + '.pkl', 'wb') as f:
            pickle.dump(obj, f)#, pickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":
    prefix = 'AtomPair_randomForest_B3LYP_supervised'
    data = loadData('model',"/Users/arindam/HOPV/")
    fpType = 'AtomPair'
    yType ='HOMO'
    funcType = 'B3LYP'
    dummySMILES = 'Cc1csc(c2ccc(c3cc(C)c(c4cc5c(s4)c(c4ccc(C)cc4)c4ccsc4c5c4ccc(C)cc4)s3)c3nsnc23)c1'
    dummyFingerprint = dummyAtomPairFromSmiles(dummySMILES)
    estimator = data['best_estimator']
    X_full = loadData("X_full_"+fpType,"/Users/arindam/HOPV/FP")
    y_true = array(loadData(yType+'_'+funcType,"/Users/arindam/HOPV/FP"))
    estimator.fit(X_full,y_true)
    saveData(estimator,'bestModel','/Users/arindam/Desktop/Website/','.')
