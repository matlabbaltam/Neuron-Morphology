import numpy as np
import os
import xml.etree.ElementTree as ET
from lxml import etree

from matplotlib import pyplot as plt

from scipy.interpolate import SmoothSphereBivariateSpline,LSQSphereBivariateSpline
#import pyvista as pv
import pymeshfix
import pickle
#import pyshtools as pysh
#import math as m
import pandas as pd
import open3d as o3d
from cell_name import sort_name,check_name
import string
import copy
import pdb
import csv



def cart2sph(coords):
    center = np.average(coords,axis = 0)
    coords_rescale = np.zeros_like(coords)
    coords_rescale[:,0] = coords[:,0]-center[0]
    coords_rescale[:,1] = coords[:,1]-center[1]
    coords_rescale[:,2] = coords[:,2]-center[2]
    x = coords_rescale[:,0]
    y = coords_rescale[:,1]
    z = coords_rescale[:,2]
    XsqPlusYsq = x**2 + y**2
    r = np.sqrt(XsqPlusYsq + z**2)               # r
    elev = np.arctan2(z,np.sqrt(XsqPlusYsq))     # [0,pi]
    elev[elev<0] = elev[elev<0]+np.pi
    az = np.arctan2(y,x)                               #[0,2pi]
    az[az<0] = az[az<0]+2*np.pi
    return r, elev, az


def sph2cart(sph):
    num_bins = len(sph[0])
    theta = np.linspace(0,np.pi,num_bins)
    phi = np.linspace(0,np.pi*2,num_bins)
    x = []
    y = []
    z = []
    for i in range(len(sph[0])):
        for j in range(len(sph[1])):
            x.append(sph[i,j]*np.sin(theta[i])*np.cos(phi[j]))
            y.append(sph[i,j]*np.sin(theta[i])*np.sin(phi[j]))
            z.append(sph[i,j]*np.cos(theta[i]))
    return x,y,z



def Xforward(dim, a, b, x, y):
    if dim == 1:
        return a[0] + x
    elif dim == 2:
        return a[0] + a[1]*x
    elif dim == 3:
        return a[0] + a[1]*x + a[2]*y
    elif dim == 4:
        return a[0] + (a[1] + a[3]*y)*x + a[2]*y
    elif dim == 5:
        return a[0] + (a[1] + a[3]*y + a[4]*x)*x + a[2]*y
    elif dim == 6:
        return a[0] + (a[1] + a[3]*y + a[4]*x)*x + (a[2] + a[5]*y)*y
    return x

def Yforward(dim, a, b, x, y):
    if dim == 1:
        return b[0] + y
    elif dim == 2:
        return b[0] + b[1]*y
    elif dim == 3:
        return b[0] + b[1]*x + b[2]*y
    elif dim == 4:
        return b[0] + (b[1] + b[3]*y)*x + b[2]*y
    elif dim == 5:
        return b[0] + (b[1] + b[3]*y + b[4]*x)*x + b[2]*y
    elif dim == 6:
        return b[0] + (b[1] + b[3]*y + b[4]*x)*x + (b[2] + b[5]*y)*y
    return y




def get_section_thick():
    file_name = "/home/tom/cell_shape/Dict-20211115T224415Z-001/Dict/sections-new.csv"
    file = open(file_name)
    csvreader = csv.reader(file)
    thick_dic = {}
    for row in csvreader:
        thick_dic[row[0]]=row[1]
    file.close()
    return thick_dic

'''
def merge_mesh(mesh):
    triangle_clusters, cluster_n_triangles, cluster_area = mesh.cluster_connected_triangles()
    triangle_clusters = np.asarray(triangle_clusters)
    cluster_n_triangles = np.asarray(cluster_n_triangles)
    cluster_area = np.asarray(cluster_area)  
    if (len(cluster_area)==1):
        return mesh
    elif:

    print (cluster_area.shape)
    pdb.set_trace()
'''
###get cell ID
def get_smooth_cell_dict(name):
    thick_dic = get_section_thick()
    file_folder = "/home/tom/cell_shape/Dict-20211115T224415Z-001/Dict/"
    #file_folder = "/home/tom/cell_shape/series-20220505T015131Z-001/series/"
    cell_name= name
    cell_names = sort_name(cell_name)
    #cell_names = ['AVG7ax']
    #print (cell_names)
    #print (cell_id.isnumeric())
    #pdb.set_trace()
    num_sections = 7671
    coords = []
    #clm = pysh.SHCoeffs.from_array
    #points = pv.Sphere().points
    #print (points)
    #pdb.set_trace()
    epsilon = 5e-10
    real_z = 0
    for serial_num in range(num_sections):
        serial_num_tmp = serial_num
        #if not os.path.exists(file_folder+'Series1.'+str(serial_num+1)+'.pkl'):
        #    serial_num_tmp = serial_num_exist
        if os.path.exists(file_folder+'Series1.'+str(serial_num_tmp+1)+'.pkl'):
            
            #serial_num_exist = serial_num_tmp
            filename = open(file_folder+'Series1.'+str(serial_num_tmp+1)+'.pkl','rb')
            tree = pickle.load(filename)
            contours = tree['contours']
            #section_prop = tree['Section']
            #thickness = float(section_prop['thickness'])
            for cell in contours:
                #if (cell_id not in cell) or ("gapjunction" in cell) or ("00syn" in cell) or ("gap junction" in cell):
                #    continue
                #mag = tree['mag']
                #if cell_id in cell:
                #    ind = cell.find(cell_id)
                #    ind_tmp = ind+len(cell_id)
                #    if len(cell)>ind_tmp and cell[ind_tmp].isdigit():
                #        continue
                for cell_id in cell_names:
                    if (cell_id.isnumeric() and cell.startswith(cell_id)) or (not cell_id.isnumeric() and cell_id in cell):
                        if not check_name(cell_id,cell):
                            continue
                        real_z = real_z+1
                        a = contours[cell]['xcoef']
                        b = contours[cell]['ycoef']
                        #if len(cell)>len(cell_id) and cell[len(cell_id)].isnumeric():
                        #    continue
                        if "R" in cell and not "R" in cell_id:
                            continue
                        if "L" in cell and not "L" in cell_id:
                            continue
                        print (serial_num,cell)
                        coordinates = contours[cell]['points']
                        coordinates = np.asarray(coordinates,dtype='float')
                        coordinates = np.reshape(coordinates,(-1,2))
                        c_p=[]
                        for i in coordinates:
                            dim = 6
                            x = i[0]
                            y = i[1]
                            u = x
                            v = y
                            x0 = 0.0           #initial guess of (x,y)
                            y0 = 0.0
                            u0 = Xforward(dim, a, b, x0,y0)          #get forward tform of initial guess
                            v0 = Yforward(dim, a, b, x0,y0)
                            i = 0
                            e = 1.0
                            while (e > epsilon) and (i < 10):
                                i += 1
                                l = a[1] + a[3]*y0 + 2.0*a[4]*x0
                                m = a[2] + a[3]*x0 + 2.0*a[5]*y0
                                n = b[1] + b[3]*y0 + 2.0*b[4]*x0
                                o = b[2] + b[3]*x0 + 2.0*b[5]*y0
                                p = l*o - m*n
                                if abs(p) > epsilon:
                                    x0 += (o*(u-u0) - m*(v-v0))/p
                                    y0 += (l*(v-v0) - n*(u-u0))/p
                                else:
                                    x0 += l*(u-u0) + n*(v-v0)
                                    y0 += m*(u-u0) + o*(v-v0)
                                u0 = Xforward(dim, a, b, x0,y0)
                                v0 = Yforward(dim, a, b, x0,y0)
                                e = abs(u-u0) + abs(v-v0)
                            result_x = x0
                            result_y = y0
                            c_p.append([result_x,result_y])
                        c_p = np.asarray(c_p,dtype='float')
                        #c_p = coordinates
                        #plt.plot(coordinates[:,0],coordinates[:,1])
                        #plt.savefig('sections/'+str(serial_num)+'.png')
                        #plt.close('all')
                        #z_dim = (real_z)*np.ones((len(coordinates),1))*float(thick_dic[str(serial_num_tmp+1)])
                        #z_dim = (serial_num+1)*np.ones((len(coordinates),1))*float(thick_dic[str(serial_num_tmp+1)])
                        #c_p = np.concatenate((c_p,z_dim),axis = 1)
                        
                        thickness = thick_dic[str(serial_num_tmp+1)]
                        num_thick = int(float(thickness)/0.01)
                        num_points = len(c_p)
                        c_p = np.repeat(c_p,num_thick,axis=0)
                        start_z = (real_z+1)*float(thick_dic[str(serial_num_tmp+1)])-float(thick_dic[str(serial_num_tmp+1)])
                        z_dim = []
                        for num in range(num_thick):
                            z_dim_tmp = list(start_z*np.ones(num_points))
                            start_z = start_z+0.01
                            z_dim = z_dim+z_dim_tmp
                        z_dim = np.reshape(np.asarray(z_dim),(-1,1))
                        
                        #print (z_dim)
                        #pdb.set_trace()
                        c_p = np.concatenate((c_p,z_dim),axis = 1)
                        #real_z = 
                        #z_dim_tmp = (real_z+1)*np.ones((len(coordinates),1))*0.01
                        #real_z = real_z+1
                        #z_dim = z_dim+z_dim_tmp
                        #z_dim = np.asarray(z_dim,dtype='float')
                        #print (z_dim.shape)
                        #pdb.set_trace()
                        #c_p = np.concatenate((c_p,z_dim),axis = 1)
                        if len(coords)==0:
                            coords = c_p
                            #print (coords)
                        else:
                            coords = np.concatenate((coords,c_p))
            #except:
                #print ('weird character in '+str(serial_num+1))
    print (coords.shape)
    #print (np.average(coords,axis=0))


    ### Interpolate in Z


    #pdb.set_trace()
    #cloud = pv.PolyData(coords)
    #print (cloud.cell_arrays)
    #surf = cloud.delaunay_2d()
    #surf.plot(show_edges=True)
    pcd_old = o3d.geometry.PointCloud()
    pcd_old.points = o3d.utility.Vector3dVector(coords)
    pcd_old.estimate_normals()
    pcd_old.orient_normals_consistent_tangent_plane(100)
    #o3d.visualization.draw_geometries([pcd_old])
    
    alpha = 2
    #radii = [0.5,0.8]
    #distances = pcd_old.compute_nearest_neighbor_distance()
    #alpha = 500 * avg_dist
    #print (alpha)
    #poisson_mesh,densities = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(pcd_old, depth=9)
    #points_2k = poisson_mesh.sample_points_poisson_disk(2000)
    poisson_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_alpha_shape(pcd_old, alpha)
    #poisson_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_ball_pivoting(pcd_old, o3d.utility.DoubleVector([radius*10,radius*20]))
    
    #triangle_clusters, cluster_n_triangles, cluster_area = poisson_mesh.cluster_connected_triangles()
    #triangle_clusters = np.asarray(triangle_clusters)
    #cluster_n_triangles = np.asarray(cluster_n_triangles)
    #cluster_area = np.asarray(cluster_area)    
    ### remove floating faces
    #num_points = len(coords) 
    #num_cluster = len(triangle_clusters)
    #print (num_cluster)

    #triangles_to_remove = cluster_n_triangles[triangle_clusters] < 0.1*num_points/num_cluster
    #poisson_mesh.remove_triangles_by_mask(triangles_to_remove)
    vertices = np.asarray(poisson_mesh.vertices)
    triangles = np.asarray(poisson_mesh.triangles)
    
    ### filling holes
    tin = pymeshfix.PyTMesh()
    tin.load_array(vertices, triangles)
    tin.fill_small_boundaries()
    vertices,triangles = tin.return_arrays()
    #vertices, triangles = pymeshfix.clean_from_arrays(vertices, triangles,verbose=False,joincomp=False,remove_smallest_components=True)
    #vertices, triangles = pymeshfix.clean_from_arrays(vertices, triangles)
    new_mesh = o3d.geometry.TriangleMesh()
    new_mesh.vertices = o3d.utility.Vector3dVector(np.asarray(vertices))
    new_mesh.triangles = o3d.utility.Vector3iVector(np.asarray(triangles).astype(np.int32))
    

    num_points = len(coords)
    triangle_clusters, cluster_n_triangles, cluster_area = new_mesh.cluster_connected_triangles()
    triangle_clusters = np.asarray(triangle_clusters)
    cluster_n_triangles = np.asarray(cluster_n_triangles)
    cluster_area = np.asarray(cluster_area)  
    #print (cluster_area.shape)
    #largest_cluster_idx = cluster_n_triangles.argsort()[-10:]
    triangles_to_remove = cluster_n_triangles[triangle_clusters] < num_points/1000
    #new_mesh.remove_triangles_by_mask(triangles_to_remove)
    #triangles_to_remove = cluster_n_triangles[triangle_clusters] < 10000
    new_mesh.remove_triangles_by_mask(triangles_to_remove)  
    #new_mesh = merge_mesh(new_mesh)
    '''
    vertices = np.asarray(new_mesh.vertices)
    triangles = np.asarray(new_mesh.triangles)
    
    
    tin = pymeshfix.PyTMesh()
    tin.load_array(vertices, triangles)
    tin.fill_small_boundaries()
    vertices,triangles = tin.return_arrays()
    meshfix = pymeshfix.MeshFix(vertices, triangles)
    meshfix.repair()
    vertices = meshfix.v
    triangles = meshfix.f
    

    new_mesh = o3d.geometry.TriangleMesh()
    new_mesh.vertices = o3d.utility.Vector3dVector(np.asarray(vertices))
    new_mesh.triangles = o3d.utility.Vector3iVector(np.asarray(triangles).astype(np.int32))
    '''
    #vertices,triangles = tin.return_arrays()
    #vertices, triangles = pymeshfix.clean_from_arrays(vertices, triangles,verbose=False,joincomp=False,remove_smallest_components=True)
    #vertices, triangles = pymeshfix.clean_from_arrays(vertices, triangles)
    #new_mesh = o3d.geometry.TriangleMesh()
    #new_mesh.vertices = o3d.utility.Vector3dVector(np.asarray(vertices))
    #new_mesh.triangles = o3d.utility.Vector3iVector(np.asarray(triangles).astype(np.int32))

    #triangle_clusters, cluster_n_triangles, cluster_area = new_mesh.cluster_connected_triangles()
    #triangle_clusters = np.asarray(triangle_clusters)
    #cluster_n_triangles = np.asarray(cluster_n_triangles)
    #cluster_area = np.asarray(cluster_area)

    ### remove floating faces
    #num_points = len(coords) 
    #num_cluster = len(triangle_clusters)
    #print (num_cluster)

    #triangles_to_remove = cluster_n_triangles[triangle_clusters] < 0.1*num_points/num_cluster
    #poisson_mesh.remove_triangles_by_mask(triangles_to_remove)

    #poisson_mesh.compute_vertex_normals()
    #poisson_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(pcd, depth=8, width=0, scale=1.1, linear_fit=False)[0]
    #poisson_mesh.compute_vertex_normals()

    return new_mesh
    

    #pcd_clean = pcd
    #cl,ind=pcd_clean.remove_radius_outlier(nb_points=10,radius=0.4)

    #o3d.visualization.draw_geometries([cl])

    #points = pv.wrap(coords)
    #print (points)
    #pv.plot(points)
    #mesh = pv.StructuredGrid(x,y,z)
    #mesh.plot(show_edges=True)


###get cell ID
def get_smooth_cell(name,alpha=1):
    file_folder = "/home/tom/cell_shape/series-20220505T015131Z-001/series/"
    cell_name= name
    cell_names = sort_name(cell_name)
    #cell_names = ['AVG7ax']
    #print (cell_names)
    #print (cell_id.isnumeric())
    #pdb.set_trace()
    num_sections = 7671
    coords = []
    #clm = pysh.SHCoeffs.from_array
    #points = pv.Sphere().points
    #print (points)
    #pdb.set_trace()
    epsilon = 5e-10
    start_z = 0
    for serial_num in range(1,num_sections):
        #if not os.path.exists(file_folder+'Series1.'+str(serial_num+1)+'.pkl'):
        #    serial_num_tmp = serial_num_exist
        if os.path.exists(file_folder+'Series1.'+str(serial_num)):
            #real_z = real_z+1
            #serial_num_exist = serial_num_tmp
            parser = etree.XMLParser(recover=True)
            with open(file_folder+'Series1.'+str(serial_num),'r',encoding="latin-1") as xml_txt:
                root = etree.fromstring((xml_txt.read()), parser = parser)
            #tree = ET.parse(file_folder+'frustration.'+str(serial_num),parser = parser)
            #root = tree.getroot()
            thickness = float(root.attrib['thickness'])
            for child in root:
                for grandchild in child:
                    if grandchild.tag=="Contour":
                        cell = grandchild.attrib['name']
                        for cell_id in cell_names:
                            if (cell_id.isnumeric() and cell.startswith(cell_id)) or (not cell_id.isnumeric() and cell_id in cell):
                                if not check_name(cell_id,cell):
                                    continue
                                #a = contours[cell]['xcoef']
                                #b = contours[cell]['ycoef']
                                #if len(cell)>len(cell_id) and cell[len(cell_id)].isnumeric():
                                #    continue
                                if "R" in cell and not "R" in cell_id:
                                    continue
                                if "L" in cell and not "L" in cell_id:
                                    continue
                                ### get transform
                                a = []
                                b = []
                                a_tmp = child.attrib['xcoef']
                                b_tmp = child.attrib['ycoef']
                                dim = int(child.attrib['dim'])
                                a_tmp = a_tmp.split(' ')
                                a_tmp = a_tmp[1:]
                                b_tmp = b_tmp.split(' ')
                                b_tmp = b_tmp[1:]
                                for a_ele in a_tmp:
                                    a.append(float(a_ele))
                                for b_ele in b_tmp:
                                    b.append(float(b_ele))

                                print (serial_num,cell)
                                if not 'points' in grandchild.attrib:
                                    continue
                                coordinates_tmp = grandchild.attrib['points']
                                coordinates_tmp = coordinates_tmp.strip()
                                coordinates_tmp = coordinates_tmp.split(',')
                                #coordinates = np.fromstring(coordinates_tmp)
                                coordinates = []
                                for ele in coordinates_tmp:
                                    ele = ele.strip()
                                    ele = ele.split(' ')
                                    if len(ele)==2:
                                        #print (ele[0],ele[1])
                                        #pdb.set_trace()
                                        coordinates.append([float(ele[0]),float(ele[1])])
                                #print (coordinates)
                                #pdb.set_trace()
                                coordinates = np.asarray(coordinates,dtype='float')
                                coordinates = np.reshape(coordinates,(-1,2))
                                #print(a,b,dim)
                                c_p=[]
                                for i in coordinates:
                                    x = i[0]
                                    y = i[1]
                                    u = x
                                    v = y
                                    x0 = 0.0           #initial guess of (x,y)
                                    y0 = 0.0
                                    u0 = Xforward(dim, a, b, x0,y0)          #get forward tform of initial guess
                                    v0 = Yforward(dim, a, b, x0,y0)
                                    i = 0
                                    e = 1.0
                                    while (e > epsilon) and (i < 10):
                                        i += 1
                                        l = a[1] + a[3]*y0 + 2.0*a[4]*x0
                                        m = a[2] + a[3]*x0 + 2.0*a[5]*y0
                                        n = b[1] + b[3]*y0 + 2.0*b[4]*x0
                                        o = b[2] + b[3]*x0 + 2.0*b[5]*y0
                                        p = l*o - m*n
                                        if abs(p) > epsilon:
                                            x0 += (o*(u-u0) - m*(v-v0))/p
                                            y0 += (l*(v-v0) - n*(u-u0))/p
                                        else:
                                            x0 += l*(u-u0) + n*(v-v0)
                                            y0 += m*(u-u0) + o*(v-v0)
                                        u0 = Xforward(dim, a, b, x0,y0)
                                        v0 = Yforward(dim, a, b, x0,y0)
                                        e = abs(u-u0) + abs(v-v0)
                                    result_x = x0
                                    result_y = y0
                                    c_p.append([result_x,result_y])
                                c_p = np.asarray(c_p,dtype='float')
                                num_thick = int(float(thickness)/0.01)
                                num_points = len(c_p)
                                c_p = np.repeat(c_p,num_thick,axis=0)
                                z_dim = []
                                for num in range(num_thick):
                                    z_dim_tmp = list(start_z*np.ones(num_points))
                                    start_z = start_z+0.01
                                    z_dim = z_dim+z_dim_tmp
                                z_dim = np.reshape(np.asarray(z_dim),(-1,1))
                                
                                #print (z_dim)
                                #pdb.set_trace()
                                c_p = np.concatenate((c_p,z_dim),axis = 1)                                
                                if len(coords)==0:
                                    coords = c_p
                                    #print (coords)
                                else:
                                    coords = np.concatenate((coords,c_p))                           
        start_z = start_z+thickness
        #except:
        #    print ('weird character in '+str(serial_num+1))
    print (coords.shape)
    #print (np.average(coords,axis=0))


    ### Interpolate in Z


    #pdb.set_trace()
    #cloud = pv.PolyData(coords)
    #print (cloud.cell_arrays)
    #surf = cloud.delaunay_2d()
    #surf.plot(show_edges=True)
    pcd_old = o3d.geometry.PointCloud()
    pcd_old.points = o3d.utility.Vector3dVector(coords)
    pcd_old.estimate_normals()
    pcd_old.orient_normals_consistent_tangent_plane(100)
    #o3d.visualization.draw_geometries([pcd_old])
    
    alpha = alpha
    #radii = [0.5,0.8]
    distances = pcd_old.compute_nearest_neighbor_distance()
    avg_dist = np.mean(distances)
    #alpha = 500 * avg_dist
    #print (alpha)
    #poisson_mesh,densities = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(pcd_old, depth=9)
    #points_2k = poisson_mesh.sample_points_poisson_disk(2000)
    poisson_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_alpha_shape(pcd_old, alpha)
    #poisson_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_ball_pivoting(pcd_old, o3d.utility.DoubleVector([radius*10,radius*20]))
    
    #triangle_clusters, cluster_n_triangles, cluster_area = poisson_mesh.cluster_connected_triangles()
    #triangle_clusters = np.asarray(triangle_clusters)
    #cluster_n_triangles = np.asarray(cluster_n_triangles)
    #cluster_area = np.asarray(cluster_area)    
    ### remove floating faces
    #num_points = len(coords) 
    #num_cluster = len(triangle_clusters)
    #print (num_cluster)

    #triangles_to_remove = cluster_n_triangles[triangle_clusters] < 0.1*num_points/num_cluster
    #poisson_mesh.remove_triangles_by_mask(triangles_to_remove)
    vertices = np.asarray(poisson_mesh.vertices)
    triangles = np.asarray(poisson_mesh.triangles)
    
    ### filling holes
    tin = pymeshfix.PyTMesh()
    tin.load_array(vertices, triangles)
    tin.fill_small_boundaries()
    vertices,triangles = tin.return_arrays()
    #vertices, triangles = pymeshfix.clean_from_arrays(vertices, triangles,verbose=False,joincomp=False,remove_smallest_components=True)
    #vertices, triangles = pymeshfix.clean_from_arrays(vertices, triangles)
    new_mesh = o3d.geometry.TriangleMesh()
    new_mesh.vertices = o3d.utility.Vector3dVector(np.asarray(vertices))
    new_mesh.triangles = o3d.utility.Vector3iVector(np.asarray(triangles).astype(np.int32))
    

    num_points = len(coords)
    triangle_clusters, cluster_n_triangles, cluster_area = new_mesh.cluster_connected_triangles()
    triangle_clusters = np.asarray(triangle_clusters)
    cluster_n_triangles = np.asarray(cluster_n_triangles)
    cluster_area = np.asarray(cluster_area)  
    #print (cluster_area.shape)
    #largest_cluster_idx = cluster_n_triangles.argsort()[-10:]
    triangles_to_remove = cluster_n_triangles[triangle_clusters] < num_points/100
    #new_mesh.remove_triangles_by_mask(triangles_to_remove)
    #triangles_to_remove = cluster_n_triangles[triangle_clusters] < 10000
    new_mesh.remove_triangles_by_mask(triangles_to_remove)
    #new_mesh.merge_close_vertices(5)  
    triangle_clusters, cluster_n_triangles, cluster_area = new_mesh.cluster_connected_triangles()
    print (len(triangle_clusters))
    print (len(cluster_n_triangles))

    if len(cluster_n_triangles)==2:
        vertices = np.asarray(new_mesh.vertices)
        triangles = np.asarray(new_mesh.triangles)
        vertices, triangles = pymeshfix.clean_from_arrays(vertices, triangles,verbose=False,joincomp=True,remove_smallest_components=False)
        new_mesh = o3d.geometry.TriangleMesh()
        new_mesh.vertices = o3d.utility.Vector3dVector(np.asarray(vertices))
        new_mesh.triangles = o3d.utility.Vector3iVector(np.asarray(triangles).astype(np.int32))

    #new_mesh = merge_mesh(new_mesh)
    
    #vertices = np.asarray(new_mesh.vertices)
    #triangles = np.asarray(new_mesh.triangles)
    
    '''
    tin = pymeshfix.PyTMesh()
    tin.load_array(vertices, triangles)
    tin.fill_small_boundaries()
    vertices,triangles = tin.return_arrays()
    meshfix = pymeshfix.MeshFix(vertices, triangles)
    meshfix.repair()
    vertices = meshfix.v
    triangles = meshfix.f
    '''
    #vertices,triangles = tin.return_arrays()
    #vertices, triangles = pymeshfix.clean_from_arrays(vertices, triangles,verbose=False,joincomp=False,remove_smallest_components=True)
    #vertices, triangles = pymeshfix.clean_from_arrays(vertices, triangles)
    #new_mesh = o3d.geometry.TriangleMesh()
    #new_mesh.vertices = o3d.utility.Vector3dVector(np.asarray(vertices))
    #new_mesh.triangles = o3d.utility.Vector3iVector(np.asarray(triangles).astype(np.int32))

    #triangle_clusters, cluster_n_triangles, cluster_area = new_mesh.cluster_connected_triangles()
    #triangle_clusters = np.asarray(triangle_clusters)
    #cluster_n_triangles = np.asarray(cluster_n_triangles)
    #cluster_area = np.asarray(cluster_area)

    ### remove floating faces
    #num_points = len(coords) 
    #num_cluster = len(triangle_clusters)
    #print (num_cluster)

    #triangles_to_remove = cluster_n_triangles[triangle_clusters] < 0.1*num_points/num_cluster
    #poisson_mesh.remove_triangles_by_mask(triangles_to_remove)

    #poisson_mesh.compute_vertex_normals()
    #poisson_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(pcd, depth=8, width=0, scale=1.1, linear_fit=False)[0]
    #poisson_mesh.compute_vertex_normals()

    return new_mesh


def get_smooth_cell_lowres(name):
    file_folder = "/home/tom/cell_shape/series_lowres-20220607T211130Z-001/series_lowres/"
    cell_name= name
    cell_names = sort_name(cell_name)
    #cell_names = ['AVG7ax']
    #print (cell_names)
    #print (cell_id.isnumeric())
    #pdb.set_trace()
    num_sections = 1868
    coords = []
    #clm = pysh.SHCoeffs.from_array
    #points = pv.Sphere().points
    #print (points)
    #pdb.set_trace()
    epsilon = 5e-10
    start_z = 0
    for serial_num in range(1,num_sections):
        #if not os.path.exists(file_folder+'Series1.'+str(serial_num+1)+'.pkl'):
        #    serial_num_tmp = serial_num_exist
        if os.path.exists(file_folder+'frustration.'+str(serial_num)):
            #real_z = real_z+1
            #serial_num_exist = serial_num_tmp
            parser = etree.XMLParser(recover=True)
            with open(file_folder+'frustration.'+str(serial_num),'r',encoding="latin-1") as xml_txt:
                root = etree.fromstring((xml_txt.read()), parser = parser)
            #tree = ET.parse(file_folder+'frustration.'+str(serial_num),parser = parser)
            #root = tree.getroot()
            thickness = float(root.attrib['thickness'])*3
            for child in root:
                for grandchild in child:
                    if grandchild.tag=="Contour":
                        cell = grandchild.attrib['name']
                        for cell_id in cell_names:
                            if (cell_id.isnumeric() and cell.startswith(cell_id)) or (not cell_id.isnumeric() and cell_id in cell):
                                if not check_name(cell_id,cell):
                                    continue
                                #a = contours[cell]['xcoef']
                                #b = contours[cell]['ycoef']
                                #if len(cell)>len(cell_id) and cell[len(cell_id)].isnumeric():
                                #    continue
                                if "R" in cell and not "R" in cell_id:
                                    continue
                                if "L" in cell and not "L" in cell_id:
                                    continue
                                ### get transform
                                a = []
                                b = []
                                a_tmp = child.attrib['xcoef']
                                b_tmp = child.attrib['ycoef']
                                dim = int(child.attrib['dim'])
                                a_tmp = a_tmp.split(' ')
                                a_tmp = a_tmp[1:]
                                b_tmp = b_tmp.split(' ')
                                b_tmp = b_tmp[1:]
                                for a_ele in a_tmp:
                                    a.append(float(a_ele))
                                for b_ele in b_tmp:
                                    b.append(float(b_ele))

                                print (serial_num,cell)
                                coordinates_tmp = grandchild.attrib['points']
                                coordinates_tmp = coordinates_tmp.strip()
                                coordinates_tmp = coordinates_tmp.split(',')
                                #coordinates = np.fromstring(coordinates_tmp)
                                coordinates = []
                                for ele in coordinates_tmp:
                                    ele = ele.strip()
                                    ele = ele.split(' ')
                                    if len(ele)==2:
                                        #print (ele[0],ele[1])
                                        #pdb.set_trace()
                                        coordinates.append([float(ele[0]),float(ele[1])])
                                #print (coordinates)
                                #pdb.set_trace()
                                coordinates = np.asarray(coordinates,dtype='float')
                                coordinates = np.reshape(coordinates,(-1,2))
                                #print(a,b,dim)
                                c_p=[]
                                for i in coordinates:
                                    x = i[0]
                                    y = i[1]
                                    u = x
                                    v = y
                                    x0 = 0.0           #initial guess of (x,y)
                                    y0 = 0.0
                                    u0 = Xforward(dim, a, b, x0,y0)          #get forward tform of initial guess
                                    v0 = Yforward(dim, a, b, x0,y0)
                                    i = 0
                                    e = 1.0
                                    while (e > epsilon) and (i < 10):
                                        i += 1
                                        l = a[1] + a[3]*y0 + 2.0*a[4]*x0
                                        m = a[2] + a[3]*x0 + 2.0*a[5]*y0
                                        n = b[1] + b[3]*y0 + 2.0*b[4]*x0
                                        o = b[2] + b[3]*x0 + 2.0*b[5]*y0
                                        p = l*o - m*n
                                        if abs(p) > epsilon:
                                            x0 += (o*(u-u0) - m*(v-v0))/p
                                            y0 += (l*(v-v0) - n*(u-u0))/p
                                        else:
                                            x0 += l*(u-u0) + n*(v-v0)
                                            y0 += m*(u-u0) + o*(v-v0)
                                        u0 = Xforward(dim, a, b, x0,y0)
                                        v0 = Yforward(dim, a, b, x0,y0)
                                        e = abs(u-u0) + abs(v-v0)
                                    result_x = x0
                                    result_y = y0
                                    c_p.append([result_x,result_y])
                                c_p = np.asarray(c_p,dtype='float')
                                num_thick = int(float(thickness)/0.01)
                                num_points = len(c_p)
                                c_p = np.repeat(c_p,num_thick,axis=0)
                                z_dim = []
                                for num in range(num_thick):
                                    z_dim_tmp = list(start_z*np.ones(num_points))
                                    start_z = start_z+0.01
                                    z_dim = z_dim+z_dim_tmp
                                z_dim = np.reshape(np.asarray(z_dim),(-1,1))
                                
                                #print (z_dim)
                                #pdb.set_trace()
                                c_p = np.concatenate((c_p,z_dim),axis = 1)                                
                                if len(coords)==0:
                                    coords = c_p
                                    #print (coords)
                                else:
                                    coords = np.concatenate((coords,c_p))                           
        start_z = start_z+thickness
        #except:
        #    print ('weird character in '+str(serial_num+1))
    print (coords.shape)
    #print (np.average(coords,axis=0))


    ### Interpolate in Z


    #pdb.set_trace()
    #cloud = pv.PolyData(coords)
    #print (cloud.cell_arrays)
    #surf = cloud.delaunay_2d()
    #surf.plot(show_edges=True)
    pcd_old = o3d.geometry.PointCloud()
    pcd_old.points = o3d.utility.Vector3dVector(coords)
    pcd_old.estimate_normals()
    pcd_old.orient_normals_consistent_tangent_plane(100)
    #o3d.visualization.draw_geometries([pcd_old])
    
    alpha = 1
    #radii = [0.5,0.8]
    distances = pcd_old.compute_nearest_neighbor_distance()
    avg_dist = np.mean(distances)
    #alpha = 500 * avg_dist
    #print (alpha)
    #poisson_mesh,densities = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(pcd_old, depth=9)
    #points_2k = poisson_mesh.sample_points_poisson_disk(2000)
    poisson_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_alpha_shape(pcd_old, alpha)
    #poisson_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_ball_pivoting(pcd_old, o3d.utility.DoubleVector([radius*10,radius*20]))
    
    #triangle_clusters, cluster_n_triangles, cluster_area = poisson_mesh.cluster_connected_triangles()
    #triangle_clusters = np.asarray(triangle_clusters)
    #cluster_n_triangles = np.asarray(cluster_n_triangles)
    #cluster_area = np.asarray(cluster_area)    
    ### remove floating faces
    #num_points = len(coords) 
    #num_cluster = len(triangle_clusters)
    #print (num_cluster)

    #triangles_to_remove = cluster_n_triangles[triangle_clusters] < 0.1*num_points/num_cluster
    #poisson_mesh.remove_triangles_by_mask(triangles_to_remove)
    vertices = np.asarray(poisson_mesh.vertices)
    triangles = np.asarray(poisson_mesh.triangles)
    
    ### filling holes
    tin = pymeshfix.PyTMesh()
    tin.load_array(vertices, triangles)
    tin.fill_small_boundaries()
    vertices,triangles = tin.return_arrays()
    #vertices, triangles = pymeshfix.clean_from_arrays(vertices, triangles,verbose=False,joincomp=False,remove_smallest_components=True)
    #vertices, triangles = pymeshfix.clean_from_arrays(vertices, triangles)
    new_mesh = o3d.geometry.TriangleMesh()
    new_mesh.vertices = o3d.utility.Vector3dVector(np.asarray(vertices))
    new_mesh.triangles = o3d.utility.Vector3iVector(np.asarray(triangles).astype(np.int32))
    

    num_points = len(coords)
    triangle_clusters, cluster_n_triangles, cluster_area = new_mesh.cluster_connected_triangles()
    triangle_clusters = np.asarray(triangle_clusters)
    cluster_n_triangles = np.asarray(cluster_n_triangles)
    cluster_area = np.asarray(cluster_area)  
    #print (cluster_area.shape)
    #largest_cluster_idx = cluster_n_triangles.argsort()[-10:]
    triangles_to_remove = cluster_n_triangles[triangle_clusters] < num_points/100
    #new_mesh.remove_triangles_by_mask(triangles_to_remove)
    #triangles_to_remove = cluster_n_triangles[triangle_clusters] < 10000
    new_mesh.remove_triangles_by_mask(triangles_to_remove)  

    #new_mesh = merge_mesh(new_mesh)
    
    #vertices = np.asarray(new_mesh.vertices)
    #triangles = np.asarray(new_mesh.triangles)
    
    '''
    tin = pymeshfix.PyTMesh()
    tin.load_array(vertices, triangles)
    tin.fill_small_boundaries()
    vertices,triangles = tin.return_arrays()
    meshfix = pymeshfix.MeshFix(vertices, triangles)
    meshfix.repair()
    vertices = meshfix.v
    triangles = meshfix.f
    '''
    #vertices,triangles = tin.return_arrays()
    #vertices, triangles = pymeshfix.clean_from_arrays(vertices, triangles,verbose=False,joincomp=False,remove_smallest_components=True)
    #vertices, triangles = pymeshfix.clean_from_arrays(vertices, triangles)
    #new_mesh = o3d.geometry.TriangleMesh()
    #new_mesh.vertices = o3d.utility.Vector3dVector(np.asarray(vertices))
    #new_mesh.triangles = o3d.utility.Vector3iVector(np.asarray(triangles).astype(np.int32))

    #triangle_clusters, cluster_n_triangles, cluster_area = new_mesh.cluster_connected_triangles()
    #triangle_clusters = np.asarray(triangle_clusters)
    #cluster_n_triangles = np.asarray(cluster_n_triangles)
    #cluster_area = np.asarray(cluster_area)

    ### remove floating faces
    #num_points = len(coords) 
    #num_cluster = len(triangle_clusters)
    #print (num_cluster)

    #triangles_to_remove = cluster_n_triangles[triangle_clusters] < 0.1*num_points/num_cluster
    #poisson_mesh.remove_triangles_by_mask(triangles_to_remove)

    #poisson_mesh.compute_vertex_normals()
    #poisson_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(pcd, depth=8, width=0, scale=1.1, linear_fit=False)[0]
    #poisson_mesh.compute_vertex_normals()

    return new_mesh

def final_mesh(cell_name):
    #smooth_surface=[]
    try:
        smooth_surface = get_smooth_cell(cell_name)
    except:
        try:
            smooth_surface = get_smooth_cell_lowres(cell_name)
        except:
            print ("Cannot find "+cell_name+" in this animal") 
    return smooth_surface 
#o3d.visualization.draw_geometries([smooth_surface_lowres],mesh_show_back_face=True)
