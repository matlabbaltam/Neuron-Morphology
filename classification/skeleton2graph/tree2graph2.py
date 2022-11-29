#import pwd
import xml.etree.ElementTree as ET
import pdb
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
#matplotlib.use('TkAgg')
from mpl_toolkits import mplot3d
import csv
import networkx as nx
from scipy.spatial import distance

### get all cell IDs
'''
IDs = []
tree = ET.parse('/home/tom/cell_shape/animaltwo/1_Ciona-ser3-part2.xml')
root = tree.getroot()
for t2_treeline in root.iter('t2_treeline'):
    print (t2_treeline.attrib["title"])
    IDs.append([t2_treeline.attrib["title"]])   
tree = ET.parse('/home/tom/cell_shape/animaltwo/1_Ciona-ser3-part-3.xml')
root = tree.getroot()
for t2_treeline in root.iter('t2_treeline'):
    IDs.append([t2_treeline.attrib["title"]])
tree = ET.parse('/home/tom/cell_shape/animaltwo/11500 Ciona Project.083018 Skeletons May4.xml')
root = tree.getroot()
for t2_treeline in root.iter('t2_treeline'):
    IDs.append([t2_treeline.attrib["title"]])
print (IDs)
fout = "Animal2_Names.csv"
with open(fout, "w") as f:
    writer = csv.writer(f)
    writer.writerows(IDs) 
pdb.set_trace()
'''

#cell_names = ["Cell 1 ACIN 2L"]
def tree2graph(cell_names):
    tree1 = ET.parse('animaltwo/1_Ciona-ser3-part2.xml')
    root1 = tree1.getroot()
    tree2 = ET.parse('animaltwo/1_Ciona-ser3-part-3.xml')
    root2 = tree2.getroot()
    tree3 = ET.parse('animaltwo/11500 Ciona Project.083018 Skeletons May4.xml')
    root3 = tree3.getroot()

    for calib in root1.iter('t2_calibration'):
        pw1 = float(calib.attrib["pixelWidth"])
        ph1 = float(calib.attrib["pixelHeight"])
        pz1 = float(calib.attrib["pixelDepth"])

    for calib in root2.iter('t2_calibration'):
        pw2 = float(calib.attrib["pixelWidth"])
        ph2 = float(calib.attrib["pixelHeight"])
        pz2 = float(calib.attrib["pixelDepth"])

    for calib in root3.iter('t2_calibration'):
        pw3 = float(calib.attrib["pixelWidth"])
        ph3 = float(calib.attrib["pixelHeight"])
        pz3 = float(calib.attrib["pixelDepth"])



    data = []

    for cell_name in cell_names:
        oid1 = {}
        for t2_layer in root1.iter('t2_layer'):
            for child in t2_layer:
                transform = child.attrib["transform"].replace('matrix','').replace('(','').replace(')','')
                transform_matrix_entries = transform.split(',')
                transform_matrix = np.zeros(9)
                i=0
                for elem in transform_matrix_entries:
                    transform_matrix[i] = float(elem)
                    i=i+1
                transform_matrix[-1]=1
                transform_matrix = np.reshape(transform_matrix,(3,3))

            #print(t2_layer.getchildren())
            #pdb.set_trace()
            oid1[t2_layer.attrib["oid"]]= [t2_layer.attrib["z"],transform_matrix]
        for t2_treeline in root1.iter('t2_treeline'):
            if t2_treeline.attrib["title"]==cell_name:
            #print (t2_treeline.attrib["title"])   
                for t2_node in t2_treeline.iter('t2_node'):
                    count = 0
                    children = t2_node.getchildren()
                    if len(children)==0:
                        print (t2_node.attrib["x"])
                    #    for child in children:
                    #        print (child.attrib["x"])
                    #for child in t2_node:
                    #    count = count+1
                    #    print (count,t2_node.attrib["x"])
                    #    if count%2==0:
                    #        pdb.set_trace()
                    #transform = t2_treeline.attrib["transform"].replace('matrix','').replace('(','').replace(')','')
                    #transform_matrix_entries = transform.split(',')
                    #transform_matrix = np.zeros(9)
                    #i=0
                    #for elem in transform_matrix_entries:
                    #    transform_matrix[i] = float(elem)
                    #    i=i+1
                    #transform_matrix[-1]=1
                    #transform_matrix = np.reshape(transform_matrix,(3,3))
                    transform_matrix = oid1[t2_node.attrib["lid"]][1]
                    x = float(t2_node.attrib["x"])
                    y = float(t2_node.attrib["y"])
                    xynew = np.matmul(transform_matrix,np.array([x,y,1]))
                    xnew = xynew[0]
                    ynew = xynew[1]
                    z = float(oid1[t2_node.attrib["lid"]][0])
                    data.append([xnew,  ynew, z ])

        oid2 = {}
        for t2_layer in root2.iter('t2_layer'):
            for child in t2_layer:
                transform = child.attrib["transform"].replace('matrix','').replace('(','').replace(')','')
                transform_matrix_entries = transform.split(',')
                transform_matrix = np.zeros(9)
                i=0
                for elem in transform_matrix_entries:
                    transform_matrix[i] = float(elem)
                    i=i+1
                transform_matrix[-1]=1
                transform_matrix = np.reshape(transform_matrix,(3,3))

            #print(t2_layer.getchildren())
            #pdb.set_trace()
            oid2[t2_layer.attrib["oid"]]= [t2_layer.attrib["z"],transform_matrix]


        for t2_treeline in root2.iter('t2_treeline'):
            #print (t2_treeline.attrib["title"])   
            for t2_node in t2_treeline.iter('t2_node'):
                if t2_treeline.attrib["title"]==cell_name:
                    transform_matrix = oid2[t2_node.attrib["lid"]][1]
                    x = float(t2_node.attrib["x"])
                    y = float(t2_node.attrib["y"])
                    xynew = np.matmul(transform_matrix,np.array([x,y,1]))
                    xnew = xynew[0]
                    ynew = xynew[1]
                    z = float(oid2[t2_node.attrib["lid"]][0])
                    data.append([xnew,  ynew, z ])


        oid3 = {}
        for t2_layer in root3.iter('t2_layer'):
            for child in t2_layer:
                transform = child.attrib["transform"].replace('matrix','').replace('(','').replace(')','')
                transform_matrix_entries = transform.split(',')
                transform_matrix = np.zeros(9)
                i=0
                for elem in transform_matrix_entries:
                    transform_matrix[i] = float(elem)
                    i=i+1
                transform_matrix[-1]=1
                transform_matrix = np.reshape(transform_matrix,(3,3))

            #print(t2_layer.getchildren())
            #pdb.set_trace()
            oid3[t2_layer.attrib["oid"]]= [t2_layer.attrib["z"],transform_matrix]
            #oid3[t2_layer.attrib["oid"]]= t2_layer.attrib["z"]
        for t2_treeline in root3.iter('t2_treeline'):
            for t2_node in t2_treeline.iter('t2_node'):
                if t2_treeline.attrib["title"]==cell_name:
                    transform_matrix = oid3[t2_node.attrib["lid"]][1]
                    x = float(t2_node.attrib["x"])
                    y = float(t2_node.attrib["y"])
                    xynew = np.matmul(transform_matrix,np.array([x,y,1]))
                    xnew = xynew[0]
                    ynew = xynew[1]
                    z = float(oid3[t2_node.attrib["lid"]][0])
                    data.append([xnew,  ynew, z ])
    data = np.asarray(data)
    print (len(data))
    G = nx.Graph()
    total_len = len(data)
    for ii in range(total_len):
        data_line = data[ii]
        coords = np.asarray(data_line)
        G.add_node(ii+1,coordinate=coords)
        if ii>0:
            point1 = data[ii-1]
            point2 = data[ii]
            weight_tmp = distance.euclidean(point1,point2)
            G.add_edge(ii-1,ii,weight=weight_tmp)
    return G

G = tree2graph(["Cell 25 MN2L"])