#!pip install scikit-bio
from skbio import DistanceMatrix
from skbio.tree import nj
import os,sys

def ReadFile(filename):
  data = open(filename, "r")
  rows = data.readlines()
  Dataset1=[]
  for i in range(0,len(rows)):
      Dataset1.append(rows[i].split())
  for i in range(0,len(Dataset1)):
    print(Dataset1[i])
  matrix=Dataset1[2:]
  ids=Dataset1[1]
  return matrix,ids

def get_par_index(text):
  s = []  
  dict_for_p = {}
  p1='('
  p2=')'
  for i, k in enumerate(text):
      if k == p1:
          s.append(i)
      if k == p2:
          dict_for_p[s.pop()] = i
  return dict_for_p

def get_clusters(text,d):
  clusters=[]
  for x,y in d.items():
      clusters.append(text[x:y+1])
  return clusters

def get_nodes(clusters):
  import re
  nodes=[]
  for i in clusters:
      i = re.findall("[a-zA-Z]+", i)
      nodes.append(i)
  return nodes

from numpy import mat
def print_nodes(matrix,ids,L,node):
    print(" ")
    R=[]
    for i in range(0,len(matrix)):
        R.append(sum(matrix[i]))

    print("Calculating Total distance from  {}".format(node))
    print(R)
    print("  ")  
    print("TD Matrix")
    for v in range(0,len(ids)):
      for d in (0,len(ids[v])):
        if ids[v]==ids[d]:
          pass
        else:
          print("({},{}): {}".format(ids[v],ids[d],matrix[v][d]),end=" ")
      print("  ")
        
    TD_matrix=[]
    for i in range(0,len(ids)):
        m=[]
        d=0
        for j in range(0,len(ids)):
            if ids[i]==ids[j]:
                adjusted_distance=0
            else:
                adjusted_distance= round(float(matrix[i][j])-((R[i]+R[j])/L-2),3)
                m.append(adjusted_distance)
        TD_matrix.append(m)
    print(" ")
    print("Calculating Net Divergence")
    for v in range(0,len(TD_matrix)):
      for d in range(0,len(TD_matrix[v])):
        if ids[v]==ids[d]:
          pass
        else:
          print("TD ({},{}): {}".format(ids[v],ids[d],TD_matrix[v][d]),end=" ")
      print("  ")
    min_distance=99999
    index1=0
    index2=0
    for i in range(0,len(TD_matrix)):
        for k in range(0,len(TD_matrix[i])):
          if TD_matrix[i][k]<min_distance:
            if ids[i]!=ids[k]:
                index1=i
                index2=k
                min_distance=TD_matrix[i][k]
     
    print(" ")
    print("lowest net Divergence",min_distance)
    print("Merging  {} and {} ".format(ids[index1],ids[index2]))

    ui= (min_distance*0.5) + (R[index1]-R[index2])/(2*L)-4
    uj= min_distance-ui
    print("Branch length of {} = {}".format(ids[index1],round(ui,2)))
    print("Branch length of {} = {}".format(ids[index2],round(uj,2)))
    print("Cluster being merged with distances between them is ")
    print("({}:{},{}:{}): {}".format(ids[index1],round(ui,2),ids[index2],round(uj,2),matrix[index1][index2]))

    new_ids=ids
    ids=[]
    TD_matrix=[]
    merged_node= new_ids[index1]+new_ids[index2]
    for i in range(0,len(new_ids)):
        if new_ids[index1]==new_ids[i] or new_ids[index2]==new_ids[i]:
          pass
        else:
          ids.append(new_ids[i])
    ids.append(merged_node)
    temp_matrix=matrix
    matrix=[]
    for i in range(0,len(new_ids)):
        lis=[]
        for m in range(0,len(new_ids)):
            i1=new_ids[i]
            i2=new_ids[m]
            if i1==i2:
                lis.append(0)
            else:
              if new_ids[i]==merged_node:
                  d1= temp_matrix[index1][m]
                  d2= temp_matrix[index2][m]
                  d3=temp_matrix[index1][index2]
                  distance=(d1+d2+d3)/2
                  lis.append(distance)
              elif new_ids[m]==merged_node:
                  d1= temp_matrix[index1][i]
                  d2= temp_matrix[index2][i]
                  d3=temp_matrix[index1][index2]
                  distance=(d1+d2+d3)/2
                  lis.append(distance)
              else:
                  first=new_ids.index(new_ids[i])
                  second=new_ids.index(new_ids[m])
                  distance=temp_matrix[first][second]
                  lis.append(distance)
        matrix.append(lis)
    print("  ")
    return matrix,ids,merged_node

def neighbor_joining(data,ids):
  new_ids=ids
  merged_node=ids[0]
  
  for i in range(0,len(data)):
      M=[float(x) for x in data[i]]
      data[i]=M
  matrix=data
  for j in range(0,len(ids)-2):
    matrix,new_ids,merged_node=print_nodes(matrix,new_ids,len(ids),merged_node)
  dm = DistanceMatrix(data, ids)
  tree = nj(dm)
  text = nj(dm, result_constructor=str)
  print("Newick format : ",text)
  par_index=get_par_index(text)
  clusters=get_clusters(text,par_index)
  nodes=get_nodes(clusters)
  print(tree.ascii_art())

file_name=sys.argv[1]
matrix,ids = ReadFile(file_name)
neighbor_joining(matrix,ids)

