import numpy as np
import networkx as nx
from scipy.spatial import distance
import pandas
import os
from matplotlib import pyplot as plt
import seaborn as sns
import pdb
from neuromorpholib.swc import load_swc, NeuronMorphology

def animal1():
    naming_file = 'name_type_bvin_split.csv'
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
            cur_cell_type = cell_type[cell_ids.index(cell_id)]

            label = cell_type_set.index(cur_cell_type)
            #print (cell_type[cell_ids==cell_id])
            labels.append(label)
            labels_str.append(cur_cell_type)
            cell_id_real.append(cell_id)
    #print (labels_str)
    all_cell_type = sorted(set(labels_str),key=labels_str.index)
    DS_edge = []
    DS_graph_indicator = []
    DS_graph_label = []
    DS_node_label= []
    DS_edge_attributes = []
    DS_node_attributes = []
    path = 'cell_skeleton/'

    num_nodes = 0
    num_edges = 0
    num_graphs = 0
    for cell in cell_id_real:
        cell_label = labels[num_graphs]
        DS_graph_label.append(str(cell_label)+'\n')
        num_graphs = num_graphs+1
        file_name = cell+'_skel_edge.obj'
        with open (path+file_name,'r') as f:
            lines  = f.readlines()
        f.close()
        edges_tmp = 0
        nodes_tmp = 0
        G = nx.Graph()
        for line in lines:
            elements = line.split(' ')
            if elements[0] == "v":
                nodes_tmp = nodes_tmp+1
                x = float(elements[1])
                y = float(elements[2])
                z = float(elements[3])
                DS_graph_indicator.append(str(num_graphs)+'\n')
                DS_node_label.append(str(1)+'\n')
                coords = np.asarray([x,y,z])
                DS_node_attributes.append(str(x)+', '+str(y)+', '+str(z)+'\n')
                G.add_node(nodes_tmp,coordinate=coords)
            if elements[0]=='l':
                node1 = int(elements[1])+num_nodes
                node2 = int(elements[2])+num_nodes
                DS_edge.append(str(node1)+', '+str(node2)+'\n')
                weight_tmp = distance.euclidean(G.nodes[node1-num_nodes]['coordinate'],G.nodes[node2-num_nodes]['coordinate'])
                DS_edge_attributes.append(str(weight_tmp)+'\n')
        num_nodes = num_nodes+nodes_tmp
    out_path = 'InfoGraph/data/neuron1/raw/'
    file1 = 'neuron1_A.txt'
    file_1_txt = open(out_path+file1,'w')
    for line in DS_edge:
        file_1_txt.write(line)
    file_1_txt.close()

    file2 = 'neuron1_graph_indicator.txt'
    file_2_txt = open(out_path+file2,'w')
    for line in DS_graph_indicator:
        file_2_txt.write(line)
    file_2_txt.close()

    file3 = 'neuron1_graph_labels.txt'
    file_3_txt = open(out_path+file3,'w')
    for line in DS_graph_label:
        file_3_txt.write(line)
    file_3_txt.close()

    file4 = 'neuron1_node_attributes.txt'
    file_4_txt = open(out_path+file4,'w')
    for line in DS_node_attributes:
        file_4_txt.write(line)
    file_4_txt.close()

    file5 = 'neuron1_node_labels.txt'
    file_5_txt = open(out_path+file5,'w')
    for line in DS_node_label:
        file_5_txt.write(line)
    file_5_txt.close()

    file6 = 'neuron1_edge_attributes.txt'
    file_6_txt = open(out_path+file6,'w')
    for line in DS_edge_attributes:
        file_6_txt.write(line)
    file_6_txt.close()
    return 0

def animal2():
    naming_file = 'animaltwo/Animal2_cell_names.csv'
    df_name = pandas.read_csv(naming_file,usecols=[0,1])
    df_name["Names of Cells"] = df_name["Names of Cells"].astype(str)
    df_name["Cell Type"] = df_name["Cell Type"].astype(str)
    cell_ids = df_name["Names of Cells"].tolist()
    cell_type = df_name["Cell Type"].tolist()
    cell_type_set = list(set(cell_type))


    path = 'swcfiles/'
    files = os.listdir(path)
    r_feat = []
    labels = []
    labels_str = []
    cell_id_real = []
    for file in files:
        cell_id = file.split(".")[0]
        cell_id = cell_id.replace('*','/')
        #xyz,r = read_skel_xyzr(path+file)
        #r_feat.append(r)
        cur_cell_type = cell_type[cell_ids.index(cell_id)]
        label = cell_type_set.index(cur_cell_type)
        #print (cell_type[cell_ids==cell_id])
        labels.append(label)
        labels_str.append(cur_cell_type)
        cell_id_real.append(cell_id)
    #print (labels_str)
    all_cell_type = sorted(set(labels_str),key=labels_str.index)
    DS_edge = []
    DS_graph_indicator = []
    DS_graph_label = []
    DS_node_label= []
    DS_edge_attributes = []
    DS_node_attributes = []

    num_nodes = 0
    num_edges = 0
    num_graphs = 0
    location = {}
    for cell in cell_id_real:
        cell_label = labels[num_graphs]
        DS_graph_label.append(str(cell_label)+'\n')
        num_graphs = num_graphs+1
        file_name = cell.replace('/','*')
        file_name = file_name+'.swc'
        with open (path+file_name,'r') as f:
            lines  = f.readlines()
        f.close()
        edges_tmp = 0
        nodes_tmp = 0
        for line in lines:
            elements = line.split(' ')
            nodes_tmp = nodes_tmp+1
            x = float(elements[2])
            y = float(elements[3])
            z = float(elements[4])
            location[elements[0]] = [x,y,z]
            DS_graph_indicator.append(str(num_graphs)+'\n')
            DS_node_label.append(str(1)+'\n')
            #coords = np.asarray([x,y,z])
            DS_node_attributes.append(str(x)+', '+str(y)+', '+str(z)+'\n')
            #G.add_node(nodes_tmp,coordinate=coords)
            node1 = int(elements[0])+num_nodes
            node2 = int(elements[6])+num_nodes
            if node2-num_nodes!=-1:
                DS_edge.append(str(node1)+', '+str(node2)+'\n')
                location_2 = location[str(node2-num_nodes)]
                weight_tmp = distance.euclidean([x,y,z],location_2)
                DS_edge_attributes.append(str(weight_tmp)+'\n')
        num_nodes = num_nodes+nodes_tmp
    out_path = 'InfoGraph/data/neuron2/raw/'
    file1 = 'neuron2_A.txt'
    file_1_txt = open(out_path+file1,'w')
    for line in DS_edge:
        file_1_txt.write(line)
    file_1_txt.close()

    file2 = 'neuron2_graph_indicator.txt'
    file_2_txt = open(out_path+file2,'w')
    for line in DS_graph_indicator:
        file_2_txt.write(line)
    file_2_txt.close()

    file3 = 'neuron2_graph_labels.txt'
    file_3_txt = open(out_path+file3,'w')
    for line in DS_graph_label:
        file_3_txt.write(line)
    file_3_txt.close()

    file4 = 'neuron2_node_attributes.txt'
    file_4_txt = open(out_path+file4,'w')
    for line in DS_node_attributes:
        file_4_txt.write(line)
    file_4_txt.close()

    file5 = 'neuron2_node_labels.txt'
    file_5_txt = open(out_path+file5,'w')
    for line in DS_node_label:
        file_5_txt.write(line)
    file_5_txt.close()

    file6 = 'neuron2_edge_attributes.txt'
    file_6_txt = open(out_path+file6,'w')
    for line in DS_edge_attributes:
        file_6_txt.write(line)
    file_6_txt.close()
    return 0


if __name__=="__main__":
    animal2()