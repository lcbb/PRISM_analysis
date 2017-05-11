# -*- coding: utf-8 -*-
"""
This script performs t-SNE on synaptic features extracted from PRISM data. 
This script loads the synaptic features from the MAT files, concatenates the 
features from the same well or from the same replicate together (depending on 
the value of "split_wells"), subsample the synapses if the sample size is too 
large, standardize the features, and gereate t-SNE maps for the features. The 
output t-SNE maps are color-coded based on the feature values. 
  
Created on Sat Nov 19 20:04:01 2016
@author: Syuan-Ming Guo
"""

import numpy as np
import os
import glob
import matplotlib.pyplot as plt
import scipy.io
from time import time
from matplotlib.ticker import NullFormatter
from sklearn import manifold
import matplotlib as mpl


mpl.rcParams['pdf.fonttype'] = 42 # change the default font settings of matplotlib
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['svg.fonttype'] = 'none'
mpl.rcParams.update({'figure.autolayout': True})
#sns.set_context("poster")
#%%
def read_features(dir_path): # load synaptic features from MAT files
    subdir_path = glob.glob(os.path.join(dir_path, '*/'))    
    subdir_name = [os.path.split(subdir[:-1])[1] for subdir in subdir_path]    
    nDir = len(subdir_name)
    file_path = os.path.join(dir_path, subdir_name[1],"mask.mat")
    mat = scipy.io.loadmat(file_path)
    IabInMarkerCell = mat['IabInMarkerCell']    
    IabInMarker = np.zeros(3,dtype=object) # define empty array for collecting features from individual MAT files
    IabWithMarker = np.zeros(3,dtype=object)
    AabWithMarker = np.zeros(3,dtype=object)
    for j in range(0, nDir): # going through each sub-folder
        file_path = os.path.join(dir_path, subdir_name[j],"mask.mat")
        mat = scipy.io.loadmat(file_path)
        IabInMarkerCell = mat['IabInMarkerCell']                 
        IabInMarkerS = np.concatenate(IabInMarkerCell[:]) # convert nested arrays to arrays 
        IabInMarkerS = np.concatenate(IabInMarkerS[:])
    
        IabWithMarkerCell = mat['IabWithMarkerCell']                
        IabWithMarkerS = np.concatenate(IabWithMarkerCell[:])
        IabWithMarkerS = np.concatenate(IabWithMarkerS[:])  
    
        AabWithMarkerCell = mat['AabWithMarkerCell']                  
        AabWithMarkerS = np.concatenate(AabWithMarkerCell[:])
        AabWithMarkerS = np.concatenate(AabWithMarkerS[:])
        
        IabInMarker[j] = IabInMarkerS # each array in "IabInMarker" contains all the features from a single well (total 3 wells)
        IabWithMarker[j] = IabWithMarkerS
        AabWithMarker[j] = AabWithMarkerS
        
    imgRoundNamesNoWO = mat['imgRoundNamesNoWO']    
    imgRoundNamesNoWO = np.concatenate(imgRoundNamesNoWO[:])
    imgRoundNamesNoWO = np.concatenate(imgRoundNamesNoWO[:])              
    return IabInMarker, IabWithMarker, AabWithMarker, imgRoundNamesNoWO

def concat_data(IabInMarkerS, IabWithMarkerS, AabWithMarkerS):
    n_features = IabInMarkerS[0].shape[1]  
    nDir = IabInMarkerS.shape[0]
    IabInMarker = np.array([]).reshape(0,n_features) # define empty array for collecting features from individual MAT files
    IabWithMarker = np.array([]).reshape(0,n_features)
    AabWithMarker = np.array([]).reshape(0,n_features)
    for j in range(0, nDir): # going through each sub-folder
        IabInMarker = np.vstack([IabInMarker, IabInMarkerS[j]])
        IabWithMarker = np.vstack([IabWithMarker, IabWithMarkerS[j]])
        AabWithMarker = np.vstack([AabWithMarker, AabWithMarkerS[j]])    
    return IabInMarker, IabWithMarker, AabWithMarker
    
def stand_feature(X):
    # standardize the feature matrix so all the features have std = 1 and min = 0, 0 is assigned to missing data  
    for j in range(0,X.shape[1]):
        Xcol = X[:,j]
        Xcol_nonzero =  Xcol[Xcol>0] ;
        out_ind = left_outlier_1d(Xcol_nonzero, -2.5) ;
        Xcol_nonzero[out_ind] = 0 ;        
        Xcol[Xcol>0] = Xcol_nonzero ;        
        Xcol = Xcol/np.std(Xcol[Xcol>0])
        Xcol[Xcol>0] = Xcol[Xcol>0] - np.amin(Xcol[Xcol>0])
        X[:,j] = Xcol
    return X
    
def is_outlier(points, thresh=3.5): # find outlier indices at both sides of the distribution
    """
    Returns a boolean array with True if points are outliers and False 
    otherwise.

    Parameters:
    -----------
        points : An numobservations by numdimensions array of observations
        thresh : The modified z-score to use as a threshold. Observations with
            a modified z-score (based on the median absolute deviation) greater
            than this value will be classified as outliers.

    Returns:
    --------
        mask : A numobservations-length boolean array.

    References:
    ----------
        Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
        Handle Outliers", The ASQC Basic References in Quality Control:
        Statistical Techniques, Edward F. Mykytka, Ph.D., Editor. 
    """
    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score > thresh
    
def left_outlier_1d(points, thresh=-3.5): # find outlier indices at only the left side of the distribution
    # find the outliers that are on the left side of the input 1D distribution only
    if len(points.shape) > 1:
        raise ValueError('left outlier is not defined for multidimentional data')    
#    points = points[:,None]
    median = np.median(points, axis=0)
    diff = points - median
    diff_abs = abs(diff)
    med_abs_deviation = np.median(diff_abs)
    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score < thresh
    
def build_features(IabInMarker, IabWithMarker, AabWithMarker,log_transform=False):# subsample the synapses if needed, stadardize the features and return the feature matrix
    IabWithMarker[:,8] = IabInMarker[:,8]
    IabWithMarker[:,12] = IabInMarker[:,12]
#
    X = np.hstack((IabWithMarker, AabWithMarker[:,0:8], AabWithMarker[:,9:12])) # leave out the areas of Tuj-1 and MAP2    
#    X = np.hstack((IabWithMarker[:,0:8], IabWithMarker[:,9:12], AabWithMarker[:,0:8], AabWithMarker[:,9:12]))
    if log_transform:
        X = np.log(X+1)    
    X = stand_feature(X)  # standardize the feature matrix
    n_sub = 10000 #subsample size
    prng = np.random.RandomState(1)
    if X.shape[0]>n_sub:
        rand_arr = prng.rand(X.shape[0])
        subsample_ind = rand_arr<(n_sub/X.shape[0])
        X = X[subsample_ind,:] # subsample the data to reduce memory usage  
#        subsample_f = math.ceil(X.shape[0]/n_sub) # (subsample fraction)^-1
#        X = X[0::subsample_f,:] # subsample the data to reduce memory usage     
    return X
    
def run_tSNE(X, n_components = 2, perplexity = 100, n_iter = 5000): 
    t0 = time()
    tsne = manifold.TSNE(n_components=n_components, perplexity=perplexity, 
                         n_iter= n_iter, method='exact',init='pca',
                         random_state=1, verbose=2)
    Y = tsne.fit_transform(X)
    t1 = time()
    print("t-SNE: %.2g sec" % (t1 - t0))
    return Y    
#%%    
def plot_tSNE(X, imgRoundNamesNoWO,perplexity, plotInd): # plot the t-SNE map
    label = imgRoundNamesNoWO # labels of targets
    labelI = ["I(" + s + ")" for s in label]
    label = np.hstack((imgRoundNamesNoWO[0:8], imgRoundNamesNoWO[9:12]))
    labelA = ["A(" + s + ")" for s in label]
    labelIA = np.hstack((labelI, labelA))
    labelIA = labelIA[plotInd]
    cmapList = ['Blues', 'Greens', 'Oranges', 'RdPu', 'Purples', 'Reds']
    
#    label = np.hstack((imgRoundNamesNoWO[0:8], imgRoundNamesNoWO[9:12]))
#    labelI = ["I(" + s + ")" for s in label]
#    labelA = ["A(" + s + ")" for s in label]
#    labelIA = np.hstack((labelI, labelA))
    #% scan colormaps
    Ind_outlier = is_outlier(Y, thresh=1.5)
#    fig = plt.figure(figsize=(20, 12))
    fig = plt.figure(figsize=(10, 6))
    plt.suptitle("t-SNE with p= %i color-coded by each feature"
                 % (perplexity), fontsize=14)
    for j in range(0, len(plotInd)):        
        color = X[~Ind_outlier,plotInd[j]]
#        ax = fig.add_subplot(4, 6, j+1)
        ax = fig.add_subplot(2, 3, j+1)
        sc = plt.scatter(Y[~Ind_outlier, 0], Y[~Ind_outlier, 1], c=color, s=1,cmap=plt.get_cmap(cmapList[j]),
                         vmin=0, vmax=6, edgecolors='none')
        plt.title("%s" % (labelIA[j])) 
        plt.colorbar(sc)
        ax.xaxis.set_major_formatter(NullFormatter())
        ax.yaxis.set_major_formatter(NullFormatter())
        ax.set_axis_bgcolor((0.7,0.7,0.7))
        plt.axis('tight')    
        plt.show()           
#%%
fig_path = 'E:/Google Drive/Python/figures/'  
plotInd = [4,6,0,1,2,3] # features that are used to color code 
for rep_num in range(1, 4): # specify replicate number    
    split_wells = False # if "True", run t-SNE for each well in the replicate, otherwise pool data from all the wells together
    
    if rep_num == 1:
        dir_path = "E:/data/Neuron/cortical/Broad_HCS/14days/20161021-crotalk in cells-t9-Rep2/post_processing_multicolor/"         
    elif rep_num == 2:        
        dir_path = "E:/data/Neuron/cortical/Broad_HCS/14days/20161027-Rep3 and Rep4/Rep3/post_processing_multicolor/"         
    elif rep_num == 3:        
        dir_path = "E:/data/Neuron/cortical/Broad_HCS/14days/20161027-Rep3 and Rep4/Rep4/post_processing_multicolor/"       
    else:
        raise ValueError('replicate number is out of range')
        
    IabInMarkerAll, IabWithMarkerAll, AabWithMarkerAll, imgRoundNamesNoWO = read_features(dir_path)
    nWell = IabInMarkerAll.shape[0]
    #%
    n_components = 2 # set the t-SNE parameters
    perplexity = 40
    n_iter = 5000
    
    if split_wells:    
        for j in range(0, nWell): # going through each well
            IabInMarker = IabInMarkerAll[j] 
            IabWithMarker = IabWithMarkerAll[j] 
            AabWithMarker = AabWithMarkerAll[j]
            X = build_features(IabInMarker, IabWithMarker, AabWithMarker)
            Y = run_tSNE(X, n_components, perplexity, n_iter)
            plot_tSNE(X, imgRoundNamesNoWO,perplexity)
            plt.savefig('E:/Google Drive/Python/figures/t-SNE_Iabwith+Aab_no_Tuj1_MAP2_p= %i_rep%i_well%i.png'% (perplexity, rep_num, j+1),dpi=300)
    #        plt.savefig('E:/Google Drive/Python/figures/t-SNE_Iabwith+Aab_no_Tuj1_MAP2_p= %i_rep%i_well%i.eps'% (perplexity, rep_num, j),dpi=300)
            plt.close("all")
    else:        
        IabInMarker, IabWithMarker, AabWithMarker = concat_data(IabInMarkerAll, IabWithMarkerAll, AabWithMarkerAll)
        X = build_features(IabInMarker, IabWithMarker, AabWithMarker, log_transform =True)
        Y = run_tSNE(X, n_components, perplexity, n_iter)
        plot_tSNE(X, imgRoundNamesNoWO,perplexity,plotInd)
#        plt.savefig('E:/Google Drive/Python/figures/t-SNE_Iabwith+Aab_new_log_exact_p= %i_rep%i.png'% (perplexity, rep_num),dpi=300)
        plt.savefig('E:/Google Drive/Python/figures/t-SNE_Iabwith+Aab_new_log_exact_sub_p= %i_rep%i.png'% (perplexity, rep_num),dpi=300)
#        plt.savefig('E:/Google Drive/Python/figures/t-SNE_Iabwith+Aab_new_log_exact_sub_p= %i_rep%i.eps'% (perplexity, rep_num),dpi=300)
        plt.savefig('E:/Google Drive/Python/figures/t-SNE_Iabwith+Aab_new_log_exact_sub_p= %i_rep%i.pdf'% (perplexity, rep_num),dpi=300)
#    plt.savefig('E:/Google Drive/Python/figures/t-SNE_Iabwith+Aab_no_Tuj1_MAP2_p= %i_rep%i.eps'% (perplexity, rep_num),dpi=300)
        plt.close("all")
#    


#%% scan perplexity values
#n_neighbors = 10
#n_components = 2
##perplexity = [10, 40, 100, 400, 1000, 4000]
##perplexity = [10, 40, 100]
#perplexity = [100]
#n_iter = [5000]
#color = X[:,0]
##%
#fig = plt.figure(figsize=(20, 8))
#plt.suptitle("t-SNE_perplexity_scan", fontsize=14)
#
#for k in range(0, len(perplexity)): 
##for k in range(0, len(n_iter)): 
#    
#
#    ax = fig.add_subplot(1, 3, k+1)
#    sc = plt.scatter(Y[:, 0], Y[:, 1], c=color, s=5,cmap=plt.cm.nipy_spectral, edgecolors='none')
#    plt.title("p= %i (%.2g sec)" % (perplexity[k],t1 - t0))
##    plt.title("n_iter= %i (%.2g sec)" % (n_iter[k],t1 - t0))
#    plt.colorbar(sc)
#    ax.xaxis.set_major_formatter(NullFormatter())
#    ax.yaxis.set_major_formatter(NullFormatter())
#    plt.axis('tight')    
#    plt.show()
##    plt.savefig('E:/Google Drive/Python/figures/t-SNE_perplexity_scan_no_Tuj1.png',dpi=300)
##    plt.savefig('E:/Google Drive/Python/figures/t-SNE_niter_scan_no_Tuj1_exact.png',dpi=300)
#    
##%% scan colormaps
#    Ind_outlier = is_outlier(Y, thresh=2)
#    fig = plt.figure(figsize=(20, 12))
#    plt.suptitle("t-SNE with p= %i color-coded by each feature"
#                 % (perplexity[k]), fontsize=14)
#    for j in range(0, X.shape[1]):        
#        color = X[~Ind_outlier,j]
#        ax = fig.add_subplot(4, 6, j+1)
#        sc = plt.scatter(Y[~Ind_outlier, 0], Y[~Ind_outlier, 1], c=color, s=1,cmap=plt.cm.nipy_spectral,
#                         vmin=0, vmax=6, edgecolors='none')
#        plt.title("%s" % (labelIA[j])) 
#        plt.colorbar(sc)
#        ax.xaxis.set_major_formatter(NullFormatter())
#        ax.yaxis.set_major_formatter(NullFormatter())
#        plt.axis('tight')    
#        plt.show()
#    plt.savefig('E:/Google Drive/Python/figures/t-SNE_color_coded_feature_Iabwith+Aab_no_Tuj1_MAP2_p= %i_rep%i.png'% (perplexity[k], rep_num),dpi=300)
#    plt.savefig('E:/Google Drive/Python/figures/t-SNE_color_coded_feature_Iabwith+Aab_no_Tuj1_MAP2_p= %i_rep%i.eps'% (perplexity[k], rep_num),dpi=300)
#    plt.savefig('E:/Google Drive/Python/figures/t-SNE_color_coded_feature_Iabwith+Aab_no_Tuj1_MAP2_p= %i_rep%i.svg'% (perplexity[k], rep_num),dpi=300)
