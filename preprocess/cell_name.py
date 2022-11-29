import pandas
import numpy as np
import math
import os
import pickle
import string

def sort_name(name):
    naming_file = 'naming.csv'
    df_name = pandas.read_csv(naming_file)
    df_name['Final published name'] = df_name['Final published name'].astype(str).str.lower()
    df_name['Alternate Name'] = df_name['Alternate Name'].astype(str).str.lower()
    df_name['Original Name'] = df_name['Original Name'].astype(str).str.lower()

    final_names = df_name['Final published name'].tolist()
    alt_names = df_name['Alternate Name'].tolist()
    orig_names = df_name['Original Name'].tolist()
    #print (final_names)
    name = name.lower()
    names=[]
    if name not in final_names:
        names.append(name)
        names.append(name.upper())
        return names
    else:
        names=[]
        for i in range(len(final_names)):
            if name == final_names[i]:
                names.append(orig_names[i].upper())
                names.append(orig_names[i])
                names.append(alt_names[i].upper())
                names.append(alt_names[i])
                names.append(name)
                names.append(name.upper())
        names = list(set(names))
        return names


def check_name(cell_name,ann_name):
    digits = string.digits
    for digit in digits:
        if (cell_name+digit) in ann_name:
            return False
        if ann_name.startswith("00"):
            return False
    return True


def name_type():
    naming_file = 'name_type.csv'
    df_name = pandas.read_csv(naming_file,usecols=[0,1])
    df_name["Cell ID"] = df_name["Cell ID"].astype(str)
    df_name["Cell Type"] = df_name["Cell Type"].astype(str)
    cell_ids = df_name["Cell ID"].tolist()
    #cell_ids.remove('nan')
    cell_type = df_name["Cell Type"].tolist()
    #cell_type.remove('nan')
    return cell_ids,cell_type

#cellids,celltype = name_type()
#print (len(cellids))
#print (len(celltype))

#name = 'pr5'
#names = sort_name(name)
#print (names)