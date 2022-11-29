from ske2graph import ske2graph, select_k
from tree2graph2 import tree2graph
import pandas
from sklearn.cluster import KMeans
import networkx as nx
import os
import numpy as np
from sklearn.ensemble import RandomForestClassifier
import pdb


def spectrum_length(spectrum):
    #spectrum = spectrum*0.0033333
    spectrum = spectrum/10**5
    if len(spectrum)==100:
        return spectrum
    else:
        if len(spectrum)>100:
            spectrum_indices = np.argpartition(spectrum,-100)[-100:]
            spectrum_truncate = spectrum[spectrum_indices]
            return spectrum_truncate
        else:
            appen_len = 100-len(spectrum)
            spectrum_all = np.concatenate((spectrum,np.zeros(appen_len)))
            return spectrum_all

def read_skel_xyzr(path):
    xyz = []
    r = []
    file = open(path,'r')
    lines = file.readlines()
    count = -1
    for line in lines:
        count = count+1
        if count>0:
            one_list = line.split(' ')
            x = float(one_list[0])
            y = float(one_list[1])
            z = float(one_list[2])
            r_one = float(one_list[3])
            xyz.append([x,y,z])
            r.append(r_one)
    return xyz,r

cell_names = ["pr19"]
Gskel2 = tree2graph(cell_names)
spectrum = nx.spectrum.laplacian_spectrum(Gskel2)
spectrum_fix = spectrum_length(spectrum)
print (spectrum_fix)
#print (spectrum_fix) 
r2_feat = [spectrum_fix]



#get all cell classes
naming_file = 'name_type.csv'
df_name = pandas.read_csv(naming_file,usecols=[0,1])
df_name["Cell ID"] = df_name["Cell ID"].astype(str)
df_name["Cell Type"] = df_name["Cell Type"].astype(str)
cell_ids = df_name["Cell ID"].tolist()
cell_type = df_name["Cell Type"].tolist()
cell_type_set = list(set(cell_type))

#print (cell_type_set)
path = 'cell_skeleton/'
files = os.listdir(path)
spec_feat = []
labels = []
labels_str = []
cell_id_real = []
for file in files:
    if "xyz_r" in file:
        cell_id = file.split("_")[0]
        filename = cell_id+'_skel_edge.obj'
        G_skel = ske2graph(filename)
        laplacian = nx.spectrum.laplacian_spectrum(G_skel)
        xyz,r = read_skel_xyzr(path+file)
        #r = laplacian
        r = np.concatenate((laplacian,np.asarray(r)))
        spec_feat.append(laplacian)
        #r_feat.append(r)
        cur_cell_type = cell_type[cell_ids.index(cell_id)]

        label = cell_type_set.index(cur_cell_type)
        #print (cell_type[cell_ids==cell_id])
        labels.append(label)
        labels_str.append(cur_cell_type)
        cell_id_real.append(cell_id)
#print (labels_str)
all_cell_type = sorted(set(labels_str),key=labels_str.index)
#print (all_cell_type)
#pdb.set_trace()

labels = []
for label_str in labels_str:
    label = all_cell_type.index(label_str)
    labels.append(label)


print (spec_feat[0])
r_feat = np.asarray(spec_feat)
print (r_feat.shape)


### cluster

kmeans = KMeans(n_clusters=80, init='k-means++',random_state=0).fit(r_feat)
print (kmeans.labels_)
cluster_labels = kmeans.labels_
cluster2label = {}
for cluster_id in cluster_labels:
    cell_indx = np.where(cluster_labels==cluster_id)[0]
    cell_indx = cell_indx.astype(int)
    #cell_labels = labels[cell_indx]
    tmp_labels = []
    for one in cell_indx:
        tmp_labels.append(labels[one])
    label = max(set(tmp_labels), key = tmp_labels.count)
    cluster2label[cluster_id] = label    
prediction_cluster = []
for cluster_id in cluster_labels:
    prediction_cluster.append(cluster2label[cluster_id])
false = 0
for i in range(len(prediction_cluster)):
    #prediction_str.append(all_cell_type[prediction[i]])
    if prediction_cluster[i]!=labels[i]:
        false=false+1
        #mis_clas.append(cell_id_real[i])
        #mis_type.append(labels_str[i])
print (1-false/len(prediction_cluster))


### RF
'''
clf = RandomForestClassifier(max_depth=5, random_state=0)
clf.fit(r_feat, labels)
prediction = clf.predict(r2_feat)
#print (r_feat[0])
print (all_cell_type[prediction[0]])
pdb.set_trace()
prediction_str = []
false = 0
mis_clas = []
mis_type = []
for i in range(len(prediction)):
    prediction_str.append(all_cell_type[prediction[i]])
    if prediction[i]!=labels[i]:
        false=false+1
        mis_clas.append(cell_id_real[i])
        mis_type.append(labels_str[i])
print (1-false/len(prediction))
print (mis_clas)
print (mis_type)
'''