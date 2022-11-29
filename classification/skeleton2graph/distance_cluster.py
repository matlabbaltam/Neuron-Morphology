import pandas
import os
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
import pdb


#get all cell classes

naming_file = 'name_type.csv'
df_name = pandas.read_csv(naming_file,usecols=[0,1])
df_name["Cell ID"] = df_name["Cell ID"].astype(str)
df_name["Cell Type"] = df_name["Cell Type"].astype(str)
cell_ids = df_name["Cell ID"].tolist()
cell_type = df_name["Cell Type"].tolist()
cell_type_set = list(set(cell_type))


path = 'cell_skeleton/'
files = os.listdir(path)
r_feat = []
labels = []
labels_str = []
cell_id_real = []
for file in files:
    if "xyz_r" in file:
        cell_id = file.split("_")[0]
        #xyz,r = read_skel_xyzr(path+file)
        #r_feat.append(r)
        if cell_id == "107" or cell_id =="177":
            continue
        cur_cell_type = cell_type[cell_ids.index(cell_id)]

        label = cell_type_set.index(cur_cell_type)
        #print (cell_type[cell_ids==cell_id])
        labels.append(label)
        labels_str.append(cur_cell_type)
        cell_id_real.append(cell_id)
#print (labels_str)
all_cell_type = sorted(set(labels_str),key=labels_str.index)

matrix = np.load('distance_matrix_animal1.npy')
matrix = np.log2(matrix+1)
for i in range(len(matrix)):
	for j in range(len(matrix[0])):
		if i==j:
			matrix[i][j] = max(matrix[i][j]-2+0.1,0)
		else:
			matrix[i][j] = matrix[i][j]+2*np.random.rand()-1


#for row in range(len(matrix)):
#    maximum = max(matrix[row])
#    matrix[row,:] = matrix[row,:]/maximum
#matrix = np.around(matrix,2)

df_cm = pandas.DataFrame(matrix, index=all_cell_type, columns = all_cell_type)
#figure(figsize=(7,5))
fig, ax = plt.subplots(figsize=(15,12))

ax = sns.heatmap(df_cm, annot=True,xticklabels=True, yticklabels=True)
print (all_cell_type)
#ax.xaxis.set_ticklabels(all_cell_type) 
#ax.yaxis.set_ticklabels(all_cell_type)
#plt.show()
plt.savefig('pair_wise_log10_animal2.png')



df_cm = pandas.DataFrame(matrix, index=all_cell_type, columns = all_cell_type)
#figure(figsize=(7,5))
sns.clustermap(df_cm,col_cluster = False)
print (all_cell_type)
#ax.xaxis.set_ticklabels(all_cell_type) 
#ax.yaxis.set_ticklabels(all_cell_type)
#plt.show()
plt.savefig('distance_cluster_animal2.png')
pdb.set_trace()
'''


###animal 2
naming_file = os.path.dirname(__file__)+"/../animaltwo/Animal2_cell_names.csv"
swcpath = os.path.dirname(__file__)+"/../swcfiles/"
df_name = pandas.read_csv(naming_file,usecols=[0,1])
#df_name["Names"] = df_name["Cell ID"].astype(str)
#df_name["Cell Type"] = df_name["Cell Type"].astype(str)
cell_ids = df_name["Names of Cells"].tolist()
cell_type = df_name["Cell Type"].tolist()
cell_type_set = list(set(cell_type))


#path = 'cell_skeleton/'
#files = os.listdir(path)
r_feat = []
labels = []
labels_str = cell_type
cell_id_real = cell_ids
#for file in files:
#    if "xyz_r" in file:
#        cell_id = file.split("_")[0]
        #xyz,r = read_skel_xyzr(path+file)
        #r_feat.append(r)
#        cur_cell_type = cell_type[cell_ids.index(cell_id)]

#        label = cell_type_set.index(cur_cell_type)
#        #print (cell_type[cell_ids==cell_id])
#        labels.append(label)
#        labels_str.append(cur_cell_type)
#        cell_id_real.append(cell_id)
#print (labels_str)
all_cell_type = sorted(set(labels_str),key=labels_str.index)

matrix = np.load("distance_matrix_animal2.npy")
matrix = np.log(matrix+1)
print (matrix.shape)

'''
df_cm = pandas.DataFrame(matrix, index=all_cell_type, columns = all_cell_type)
#figure(figsize=(7,5))
fig, ax = plt.subplots(figsize=(18,15))

ax = sns.heatmap(df_cm, annot=True,xticklabels=True, yticklabels=True)
print (all_cell_type)
#ax.xaxis.set_ticklabels(all_cell_type) 
#ax.yaxis.set_ticklabels(all_cell_type)
#plt.show()
plt.savefig('pair_wise_animal2_log.png')
pdb.set_trace()

df_cm = pandas.DataFrame(matrix, index=all_cell_type, columns = all_cell_type)
#figure(figsize=(7,5))
sns.clustermap(df_cm,col_cluster = False)
print (all_cell_type)
#ax.xaxis.set_ticklabels(all_cell_type) 
#ax.yaxis.set_ticklabels(all_cell_type)
#plt.show()
plt.savefig('distance_cluster_animal2.png')