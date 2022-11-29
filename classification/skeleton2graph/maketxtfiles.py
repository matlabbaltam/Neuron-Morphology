import xml.etree.ElementTree as ET
import pdb
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
#matplotlib.use('TkAgg')
from mpl_toolkits import mplot3d
#import navis
import pandas
from neuromorpholib.swc import load_swc, NeuronMorphology
import os

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


def tree2swc(cell_names):
    tree1 = ET.parse('animaltwo/1_Ciona-ser3-part2 May31 inversion.xml')
    root1 = tree1.getroot()
    tree2 = ET.parse('animaltwo/1_Ciona-ser3-part-3 inversion May 31.xml')
    root2 = tree2.getroot()
    tree3 = ET.parse('animaltwo/11500 Ciona Project.083018 May 31.xml')
    root3 = tree3.getroot()

    for calib in root1.iter('t2_calibration'):
        pw1 = float(calib.attrib["pixelWidth"])
        ph1 = float(calib.attrib["pixelHeight"])

    for calib in root2.iter('t2_calibration'):
        pw2 = float(calib.attrib["pixelWidth"])
        ph2 = float(calib.attrib["pixelHeight"])

    for calib in root3.iter('t2_calibration'):
        pw3 = float(calib.attrib["pixelWidth"])
        ph3 = float(calib.attrib["pixelHeight"])



    data = []

    for cell_name in cell_names:
        swc_list = []
        nodeID = 1

        oid1 = {}
        for t2_layer in root1.iter('t2_layer'):
            oid1[t2_layer.attrib["oid"]]= t2_layer.attrib["z"]
        
        for t2_treeline in root1.iter('t2_treeline'):
            if t2_treeline.attrib["title"]==cell_name:
                for t2_node in t2_treeline.iter('t2_node'):
                    t2_node.attrib['__ID__'] = nodeID
                    children = t2_node.getchildren()
                    for child in children:
                        child.attrib['__parent__'] = nodeID
                    nodeID = nodeID+1
        for t2_treeline in root1.iter('t2_treeline'):
            if t2_treeline.attrib["title"]==cell_name:
                ### affine transform for x,y
                transform = t2_treeline.attrib["transform"].replace('matrix','').replace('(','').replace(')','')
                transform_matrix_entries = transform.split(',')
                transform_matrix = np.zeros((3,3))
                transform_matrix[0,0] = float(transform_matrix_entries[0])
                transform_matrix[1,0] = float(transform_matrix_entries[1])
                transform_matrix[0,1] = float(transform_matrix_entries[2])
                transform_matrix[1,1] = float(transform_matrix_entries[3])
                transform_matrix[0,2] = float(transform_matrix_entries[4])
                transform_matrix[1,2] = float(transform_matrix_entries[5])
                transform_matrix[2,2] = 1
                for t2_node in t2_treeline.iter('t2_node'):
                    swc_line = ''
                    swc_line = swc_line+str(t2_node.attrib['__ID__'])+' '+'0'+' '
                    x = float(t2_node.attrib["x"])
                    y = float(t2_node.attrib["y"])
                    xynew = np.matmul(transform_matrix,np.array([x,y,1]))
                    xnew = xynew[0]*pw1
                    ynew = xynew[1]*ph1
                    z = float(oid1[t2_node.attrib["lid"]])*pw1
                    swc_line = swc_line+str(xnew)+ ' '+str(ynew)+' '+str(z)+' '+'0'+' '
                    if '__parent__' not in t2_node.attrib:
                        swc_line = swc_line+'-1'+'\n'
                    else:
                        swc_line = swc_line+str(t2_node.attrib['__parent__'])+'\n'
                    swc_list.append(swc_line)
                    #print (transform_matrix)
                    
        
        oid2 = {}
        for t2_layer in root2.iter('t2_layer'):
            oid2[t2_layer.attrib["oid"]]= t2_layer.attrib["z"]
        
        for t2_treeline in root2.iter('t2_treeline'):
            if t2_treeline.attrib["title"]==cell_name:
                for t2_node in t2_treeline.iter('t2_node'):
                    t2_node.attrib['__ID__'] = nodeID
                    children = t2_node.getchildren()
                    for child in children:
                        child.attrib['__parent__'] = nodeID
                    nodeID = nodeID+1
        for t2_treeline in root2.iter('t2_treeline'):
            if t2_treeline.attrib["title"]==cell_name:
                ### affine transform for x,y
                transform = t2_treeline.attrib["transform"].replace('matrix','').replace('(','').replace(')','')
                transform_matrix_entries = transform.split(',')
                transform_matrix = np.zeros((3,3))
                transform_matrix[0,0] = float(transform_matrix_entries[0])
                transform_matrix[1,0] = float(transform_matrix_entries[1])
                transform_matrix[0,1] = float(transform_matrix_entries[2])
                transform_matrix[1,1] = float(transform_matrix_entries[3])
                transform_matrix[0,2] = float(transform_matrix_entries[4])
                transform_matrix[1,2] = float(transform_matrix_entries[5])
                transform_matrix[2,2] = 1
                for t2_node in t2_treeline.iter('t2_node'):
                    swc_line = ''
                    swc_line = swc_line+str(t2_node.attrib['__ID__'])+' '+'0'+' '
                    x = float(t2_node.attrib["x"])
                    y = float(t2_node.attrib["y"])
                    xynew = np.matmul(transform_matrix,np.array([x,y,1]))
                    xnew = xynew[0]*pw2
                    ynew = xynew[1]*ph2
                    z = float(oid2[t2_node.attrib["lid"]])*pw2
                    swc_line = swc_line+str(xnew)+ ' '+str(ynew)+' '+str(z)+' '+'0'+' '
                    if '__parent__' not in t2_node.attrib:
                        swc_line = swc_line+'-1'+'\n'
                    else:
                        swc_line = swc_line+str(t2_node.attrib['__parent__'])+'\n'
                    swc_list.append(swc_line)

        oid3 = {}
        for t2_layer in root3.iter('t2_layer'):
            oid3[t2_layer.attrib["oid"]]= t2_layer.attrib["z"]
        
        for t2_treeline in root3.iter('t2_treeline'):
            if t2_treeline.attrib["title"]==cell_name:
                for t2_node in t2_treeline.iter('t2_node'):
                    t2_node.attrib['__ID__'] = nodeID
                    children = t2_node.getchildren()
                    for child in children:
                        child.attrib['__parent__'] = nodeID
                    nodeID = nodeID+1
        for t2_treeline in root3.iter('t2_treeline'):
            if t2_treeline.attrib["title"]==cell_name:
                ### affine transform for x,y
                transform = t2_treeline.attrib["transform"].replace('matrix','').replace('(','').replace(')','')
                transform_matrix_entries = transform.split(',')
                transform_matrix = np.zeros((3,3))
                transform_matrix[0,0] = float(transform_matrix_entries[0])
                transform_matrix[1,0] = float(transform_matrix_entries[1])
                transform_matrix[0,1] = float(transform_matrix_entries[2])
                transform_matrix[1,1] = float(transform_matrix_entries[3])
                transform_matrix[0,2] = float(transform_matrix_entries[4])
                transform_matrix[1,2] = float(transform_matrix_entries[5])
                transform_matrix[2,2] = 1
                for t2_node in t2_treeline.iter('t2_node'):
                    swc_line = ''
                    swc_line = swc_line+str(t2_node.attrib['__ID__'])+' '+'0'+' '
                    x = float(t2_node.attrib["x"])
                    y = float(t2_node.attrib["y"])
                    xynew = np.matmul(transform_matrix,np.array([x,y,1]))
                    xnew = xynew[0]*pw3
                    ynew = xynew[1]*ph3
                    z = float(oid3[t2_node.attrib["lid"]])*pw3
                    swc_line = swc_line+str(xnew)+ ' '+str(ynew)+' '+str(z)+' '+'0'+' '
                    if '__parent__' not in t2_node.attrib:
                        swc_line = swc_line+'-1'+'\n'
                    else:
                        swc_line = swc_line+str(t2_node.attrib['__parent__'])+'\n'
                    swc_list.append(swc_line)
    if len(swc_list)==0:
        return cell_names
    path = "swcfiles/"
    filename = '_'.join(cell_names)+'.swc'
    fileswc = open(path+filename,'w')
    for line in swc_list:
        fileswc.write(line)
    fileswc.close()
    '''
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
    '''
    return 0


