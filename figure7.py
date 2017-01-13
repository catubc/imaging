#from tsf_ptcs_classes import *
#from distributions import *
#from sequential_firing import *
#from sta_utils import *

from scipy import stats
from scipy import signal

import numpy as np
import time, math
import sys
import os.path
import multiprocessing as mp
import re

from pylab import *
from scipy.interpolate import interp1d
import struct, array, csv
import scipy.optimize 
from scipy.optimize import curve_fit
import pandas as pd
from scipy.interpolate import UnivariateSpline
import matplotlib.mlab as mlab
from matplotlib import animation
from PIL import Image


import pylab as m
import mpl_toolkits.mplot3d as m3

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import pandas #comment this out or move one line down - and it works!

main_dir = '/media/cat/8TB/in_vivo/tim/dongsheng/clusters/'

#import mdp
colors=['black','blue','indianred','lightseagreen','lightsalmon','dodgerblue','mediumvioletred','indianred','lightsalmon','pink','darkolivegreen']
#colors=['black','blue','lightseagreen','lightsalmon','dodgerblue','mediumvioletred','indianred','lightsalmon','pink','darkolivegreen']


def PCA(X, n_components):
    from sklearn import decomposition

    pca = decomposition.PCA(n_components)
    pca.fit(X)
    X=pca.transform(X)

    xx=np.zeros(len(X))
    yy=np.zeros(len(X))
    zz=np.zeros(len(X))

    coords = []
    for i in range(len(X)):
         coords.append([X[i][0], X[i][1], X[i][2]])

    pca_data = np.zeros((len(X),3),dtype=np.float32)
    for i in range(len(X)):
        pca_data[i][0]=X[i][0]
        pca_data[i][1]=X[i][1]
        pca_data[i][2]=X[i][2]
    
    return pca_data, np.array(coords).T

temp = []
color_index=[]

cluster_ids_cortex = [0,1,2]
location = []
cell_indexes=[]

for c in cluster_ids_cortex:
    temp.append(np.loadtxt(main_dir + 'cluster_state_0_class_'+str(c)))
    n_traces = len(np.loadtxt(main_dir +'cluster_state_0_class_'+str(c)))
    cell_indexes.extend(np.loadtxt(main_dir +'cluster_state_0_class_'+str(c)+"_index"))
    for l in range(n_traces):
        location.append('cortex')

print len(temp[0]),len(temp[1]),len(temp[2])

cluster_ids_cortex_subcortex = [0,1,2,3,4]
for c in cluster_ids_cortex_subcortex:
    temp.append(np.loadtxt(main_dir +'cluster_state_2_class_'+str(c)))
    n_traces = len(np.loadtxt(main_dir +'cluster_state_2_class_'+str(c)))
    cell_indexes.extend(np.loadtxt(main_dir +'cluster_state_2_class_'+str(c)+"_index"))

    for l in range(n_traces):
        location.append('subcortical')

clusters =[]
for c in range(len(temp)):
    for p in range(len(temp[c])):
        clusters.append(temp[c][p])
        color_index.append(colors[c])

cell_indexes=np.array(cell_indexes)
print cell_indexes

clusters = np.array(clusters)
print clusters.shape
print cell_indexes.shape
print len(location)
#quit()

#clusters = np.array(clusters)
#*********** PCA AND CLUSTERING OF DATA *********************
#PCA

X = clusters
n_components = 3
pca_data, coords = PCA(X, n_components)

#********* KMEANS CLUSTERING *******
if True:
    n_clusters_cortex = 3
    #cluster_labels = KMEANS(pca_data, n_clusters_cortex)

    from sklearn import cluster, datasets
    clustering = cluster.KMeans(n_clusters_cortex, max_iter=1000, random_state = 1032)
    clustering.fit(pca_data)
    cluster_labels = clustering.labels_
    
    cell_kmeans_indexes=[]
    color_index = []
    for i in range(len(X)):
        color_index.append(colors[cluster_labels[i]])
        cell_kmeans_indexes.append(cluster_labels[i])

print cell_kmeans_indexes
print len(cell_kmeans_indexes)

for n in range(n_clusters_cortex):
    class_pca_cortex = []
    class_pca_cortex_index = []
    class_pca_subcortical = []  
    class_pca_subcortical_index = []
    for i in range(len(cluster_labels)):
        if cluster_labels[i]==n: #If cluster label value matches current class 
            if location[i]=='cortex':
                class_pca_cortex.append(clusters[i])
                class_pca_cortex_index.append(cell_indexes[i])  #Save original index of trace 
            else:
                class_pca_subcortical.append(clusters[i])
                class_pca_subcortical_index.append(cell_indexes[i])
    
    np.savetxt(main_dir +'classout_'+str(n)+'_state_0',class_pca_cortex)
    np.savetxt(main_dir +'classout_'+str(n)+'_state_0_index',class_pca_cortex_index)

    np.savetxt(main_dir +'classout_'+str(n)+'_state_2',class_pca_subcortical)
    np.savetxt(main_dir +'classout_'+str(n)+'_state_2_index',class_pca_subcortical_index)


#PLOT PCA CLUSTERS
cmhot = plt.get_cmap("hot")
#fig = plt.figure()
#ax = fig.gca(projection='3d')
#ax.scatter(coords[0][:224],coords[1][:224],coords[2][:224], s=120, c=color_index[:224], edgecolor=color_index[:224],marker='s')
#ax.scatter(coords[0][224:],coords[1][224:],coords[2][224:], s=120, c=color_index[224:], edgecolor=color_index[224:],marker=">")

ax = plt.subplot(111)
ax.scatter(coords[0][:224],coords[1][:224], s=120, c=color_index[:224], edgecolor=color_index[:224],marker='o', linewidth='2', alpha=0.8)
ax.scatter(coords[0][224:],coords[1][224:], s=120, c='none', edgecolor=color_index[224:],marker="o", linewidth='2', alpha=0.8)


#ax.set_title("Awake Cortex, patterns: "+str(cluster_ids[0])+" and "+str(cluster_ids[1]))


plt.show()

if True:
    bar_cortex = []
    bar_subcortical = []
    for c in range(n_clusters_cortex):
        ax=plt.subplot(2,3,1+c)
        ave_shape = np.zeros(len(clusters[0]) , dtype=np.float32)
        traces = []
        cell_counter=0
        cortex_counter = 0
        for i in range(len(cluster_labels)):
            if cluster_labels[i]==c:
                #plt.plot(clusters[i], color='black', alpha=0.1)
                vmax = np.max(np.abs(clusters[i]))
                traces.append(clusters[i]/vmax)
                cell_counter+=1
                if location[i]=='cortex': cortex_counter+=1

        trace_ave = np.average(traces, axis=0)
        trace_std = np.std(traces, axis=0)
        plt.plot(trace_ave, color=colors[c], linewidth=5, alpha=1)
        t = np.arange(0, len(trace_ave), 1)
        ax.fill_between(t, trace_ave+trace_std, trace_ave-trace_std, color=colors[c], alpha=0.25)

        plt.plot([300,300],[-1,1], color='black', linewidth =1)
        plt.plot([0,600],[0,0], color='black', linewidth =1)
        
        plt.title("Cortex: "+ str(cortex_counter)+ " Subcortex: "+ str(cell_counter-cortex_counter)) 
        
        plt.ylim(-1,1)
        
        bar_cortex.append(cortex_counter)
        bar_subcortical.append(cell_counter-cortex_counter)

    plt.show()

#***************************************************
import matplotlib.patches as mpatches

bar_cortex=[]
bar_subcortical = []
for n in range(3):
    temp = np.int32(np.loadtxt(main_dir +'classout_'+str(n)+'_state_0_index'))
    bar_cortex.append(len(temp))

for n in range(3):
    temp = np.int32(np.loadtxt(main_dir +'classout_'+str(n)+'_state_2_index'))
    bar_subcortical.append(len(temp))

bar_cortex[0], bar_cortex[1] = bar_cortex[1], bar_cortex[0]
bar_subcortical[0], bar_subcortical[1] = bar_subcortical[1], bar_subcortical[0]
colors[0], colors[1] = colors[1], colors[0]

xx = np.arange(n_clusters_cortex)

print len(xx), len(bar_cortex)
ax=plt.subplot(111)
plt.bar(xx, bar_cortex, .9, color=colors[:n_clusters_cortex], linewidth=2, alpha=.5)

hatch = ['///', '*', '---']
for j in range(1, 2,1):
    plt.bar(xx, bar_subcortical, .9, color=colors[:n_clusters_cortex], bottom = bar_cortex, hatch='/', linewidth=2, alpha=.5)

#Plot pie-charts
solid = mpatches.Patch(facecolor = 'white', edgecolor="black")
hashed = mpatches.Patch(facecolor='white', edgecolor="black", hatch='//')

labels = ['Cortex', 'Subcortical']

ax.legend([solid, hashed], labels, fontsize=12, loc=0)

plt.ylabel("# STMTDs", fontsize=30)
ax.tick_params(axis='both', labelsize=30)
plt.show()



