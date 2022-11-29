import numpy as np
import networkx as nx
from scipy.spatial import distance
import pandas
import os
from matplotlib import pyplot as plt
import seaborn as sns
import pdb
from neuromorpholib.swc import load_swc, NeuronMorphology

def ske2graph(filename):
    path = 'cell_skeleton/'
    #filename = '1_skel_edge.obj'
    G = nx.Graph()
    with open (path+filename,'r') as f:
        lines  = f.readlines()
    f.close()
    point = 0
    for line in lines:
        elements = line.split(' ')
        if elements[0] == "v":
            point = point+1
            x = float(elements[1])
            y = float(elements[2])
            z = float(elements[3])
            coords = np.asarray([x,y,z])
            G.add_node(point,coordinate=coords)
        if elements[0]=='l':
            node1 = int(elements[1])
            node2 = int(elements[2])
            weight_tmp = distance.euclidean(G.nodes[node1]['coordinate'],G.nodes[node2]['coordinate'])
            G.add_edge(node1,node2,weight = weight_tmp)
    return G

def select_k(spectrum, minimum_energy = 0.9):
    running_total = 0.0
    total = sum(spectrum)
    if total == 0.0:
        return len(spectrum)
    for i in range(len(spectrum)):
        running_total += spectrum[i]
        if running_total / total >= minimum_energy:
            return i + 1
    return len(spectrum)

def jaccard_similarity(g, h):
    i = set(g).intersection(h)
    return round(len(i) / (len(g) + len(h) - len(i)),3)


def pairwise_cluster_animal1():
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
    #print (all_cell_type)
    #pdb.set_trace()
    cell_type_dict = dict((x,[]) for x in all_cell_type)
    #print (cell_type_dict.items())

    for i in range(len(cell_id_real)):
        cell_id = cell_id_real[i]
        cell_type = labels_str[i]
        #cell_type_dict[cell_type] = cell_type_dict[cell_type] +[cell_id]
        cell_type_dict[cell_type].append(cell_id)
    print (cell_type_dict['MGIN'])
    print (len(cell_type_dict))
    pdb.set_trace()


    #labels = []
    #for label_str in labels_str:
    #    label = all_cell_type.index(label_str)
    #    labels.append(label)

    matrix = np.zeros((len(all_cell_type),len(all_cell_type)))
    for cell_type_1 in all_cell_type:
        cell_1_list = cell_type_dict[cell_type_1]
        for cell_type_2 in all_cell_type:
            cell_2_list = cell_type_dict[cell_type_2]

            ### between class cell shape distance
            if cell_type_1!=cell_type_2:
                sim_tmp = 0
                for cell_1 in cell_1_list:
                    for cell_2 in cell_2_list:
                        file_1 = cell_1+'_skel_edge.obj'
                        file_2 = cell_2+'_skel_edge.obj'
                        G1 = ske2graph(file_1)
                        G2 = ske2graph(file_2)
                        laplacian1 = nx.spectrum.laplacian_spectrum(G1)
                        laplacian2 = nx.spectrum.laplacian_spectrum(G2)
                        k1 = select_k(laplacian1)
                        k2 = select_k(laplacian2)
                        k = max(k1, k2)
                        similarity = sum((laplacian1[:k] - laplacian2[:k])**2)
                        #print (cell_1+'and'+cell_2+':'+str(similarity))
                        sim_tmp = sim_tmp+similarity
                sim_tmp = sim_tmp/len(cell_1_list)/len(cell_2_list)
                rows = all_cell_type.index(cell_type_1)
                columns = all_cell_type.index(cell_type_2)
                matrix[rows,columns] = sim_tmp
                print (cell_type_1+' and '+cell_type_2+' average: '+str(sim_tmp))
            else:
                sim_tmp = 0
                ii = len(cell_1_list)
                print (len(cell_2_list))
                if ii==1:
                    continue
                for i in range(ii-1):
                    for j in range(i+1,ii):
                        cell_1 = cell_1_list[i]
                        cell_2 = cell_2_list[j]
                        file_1 = cell_1+'_skel_edge.obj'
                        file_2 = cell_2+'_skel_edge.obj'
                        G1 = ske2graph(file_1)
                        G2 = ske2graph(file_2)
                        laplacian1 = nx.spectrum.laplacian_spectrum(G1)
                        laplacian2 = nx.spectrum.laplacian_spectrum(G2)
                        k1 = select_k(laplacian1)
                        k2 = select_k(laplacian2)
                        k = min(k1, k2)
                        similarity = sum((laplacian1[:k] - laplacian2[:k])**2)
                        sim_tmp = sim_tmp+similarity
                        #print (cell_1+'and'+cell_2+':'+str(similarity))
                sim_tmp = sim_tmp/len(cell_1_list)/(len(cell_2_list)-1)
                rows = all_cell_type.index(cell_type_1)
                columns = all_cell_type.index(cell_type_2)
                matrix[rows,columns] = sim_tmp
                print (cell_type_1+' and ' +cell_type_2+' average: '+str(sim_tmp))
    print (matrix.shape)

    np.save("distance_matrix_animal1.npy",matrix)

def pairwise_cluster_animal2():
    #get all cell classes
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
    #print (all_cell_type)
    #pdb.set_trace()
    cell_type_dict = dict((x,[]) for x in all_cell_type)
    #print (cell_type_dict.items())

    for i in range(len(cell_id_real)):
        cell_id = cell_id_real[i]
        cell_type = labels_str[i]
        #cell_type_dict[cell_type] = cell_type_dict[cell_type] +[cell_id]
        cell_type_dict[cell_type].append(cell_id)
    print (cell_type_dict['MGIN'])
    

    #labels = []
    #for label_str in labels_str:
    #    label = all_cell_type.index(label_str)
    #    labels.append(label)

    matrix = np.zeros((len(all_cell_type),len(all_cell_type)))
    for cell_type_1 in all_cell_type:
        cell_1_list = cell_type_dict[cell_type_1]
        for cell_type_2 in all_cell_type:
            cell_2_list = cell_type_dict[cell_type_2]

            ### between class cell shape distance
            if cell_type_1!=cell_type_2:
                sim_tmp = 0
                for cell_1 in cell_1_list:
                    for cell_2 in cell_2_list:
                        cell_1 = cell_1.replace('/','-')
                        #cell_2_split = cell_2.split('/')
                        cell_2 = cell_2.replace('/','-')
                        print (cell_1+ ' and ' +cell_2)
                        file_1 = swcpath+ cell_1+'.swc'
                        file_2 = swcpath+ cell_2+'.swc'
                        if (not os.path.exists(file_1)) or (not os.path.exists(file_2)):
                            continue
                        G1 = load_swc(file_1).get_graph().to_undirected()
                        G2 = load_swc(file_2).get_graph().to_undirected()
                        laplacian1 = nx.spectrum.laplacian_spectrum(G1)
                        laplacian2 = nx.spectrum.laplacian_spectrum(G2)
                        k1 = select_k(laplacian1)
                        k2 = select_k(laplacian2)
                        k = min(k1, k2)
                        similarity = sum((laplacian1[:k] - laplacian2[:k])**2)
                        print (cell_1+' and '+cell_2+':'+str(similarity))
                        sim_tmp = sim_tmp+similarity
                sim_tmp = sim_tmp/len(cell_1_list)/len(cell_2_list)
                rows = all_cell_type.index(cell_type_1)
                columns = all_cell_type.index(cell_type_2)
                matrix[rows,columns] = sim_tmp
                print (cell_type_1+' and '+cell_type_2+' average: '+str(sim_tmp))
            else:
                sim_tmp = 0
                ii = len(cell_1_list)
                print (len(cell_2_list))
                if ii==1:
                    continue
                for i in range(ii-1):
                    for j in range(i+1,ii):
                        cell_1 = cell_1_list[i]
                        cell_2 = cell_2_list[j]
                        cell_1 = cell_1.replace('/','-')
                        cell_2 = cell_2.replace('/','-')
                        print (cell_1+' and '+cell_2)
                        file_1 = swcpath+ cell_1+'.swc'
                        file_2 = swcpath+ cell_2+'.swc'
                        if (not os.path.exists(file_1)) or (not os.path.exists(file_2)):
                            continue
                        G1 = load_swc(file_1).get_graph().to_undirected()
                        G2 = load_swc(file_2).get_graph().to_undirected()
                        laplacian1 = nx.spectrum.laplacian_spectrum(G1)
                        laplacian2 = nx.spectrum.laplacian_spectrum(G2)
                        k1 = select_k(laplacian1)
                        k2 = select_k(laplacian2)
                        k = min(k1, k2)
                        similarity = sum((laplacian1[:k] - laplacian2[:k])**2)
                        sim_tmp = sim_tmp+similarity
                        print (cell_1+' and '+cell_2+':'+str(similarity))
                sim_tmp = sim_tmp/len(cell_1_list)/(len(cell_2_list)-1)
                rows = all_cell_type.index(cell_type_1)
                columns = all_cell_type.index(cell_type_2)
                matrix[rows,columns] = sim_tmp
                print (cell_type_1+' and ' +cell_type_2+' average: '+str(sim_tmp))
    print (matrix)

    np.save("distance_matrix_animal2.npy",matrix)


pairwise_cluster_animal1()