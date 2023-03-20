
import biotite.sequence.phylo as phylo
import biotite.sequence.graphics as graphics
import numpy as np
import matplotlib.pyplot as plt
import os,sys
import re

def ReadFile(filename):
  data = open(filename, "r")
  rows = data.readlines()
  Dataset1=[]
  for i in range(0,len(rows)):
      Dataset1.append(rows[i].split())
  matrix=Dataset1[2:]
  ids=Dataset1[1]
  return matrix,ids

def char_replacement(text, dic):
    for x, y in dic.items():
        text = text.replace(str(x), y)
    return text

def get_parentheses_index(text):
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

def replace_index(ids,d,text):
  clusters=[]
  for x,y in d.items():
      clusters.append(text[x:y+1])
  length=len(ids)
  replace_index={}
  for i in range(0,length):
      replace_index[i]=ids[i]
  for i in range(0,len(clusters)):
      clusters[i]=char_replacement(clusters[i],replace_index)
  return clusters,replace_index

def display_cluster(clusters):
  print("Valid Distance Matrix Loaded")
  for i in range(0,len(clusters)):
    print("cluster {} : {} = {}".format(i+1,clusters[i].replace(",","+"),clusters[i]))
  print("**************************FINAL P-TREE**************************")

def UPGMA(dataset,ids):
  data=np.array(dataset).astype(int)
  tree = phylo.upgma(data)
  text=tree.to_newick(include_distance=False)
  dict_for_index=get_parentheses_index(text)
  clusters,dic = replace_index(ids,dict_for_index,text)

  n=tree.to_newick(include_distance=False)
  for i in range(0,len(ids)):
      k=str(i)
      n=n.replace(k,ids[i])
  print('newick Format : ',n)
    
  newick_format=tree.to_newick(include_distance=True)
  for i in range(0,len(ids)):
      k=str(i)+":"
      l=ids[i]+":"
      newick_format=newick_format.replace(k,l)
  print("Tree with Branch lengths : ",newick_format)
  distances= (re.findall("\d+\.\d+",newick_format))
  res = [float(ele) for ele in distances]
  print("Final Average distance : ",sum(res)/len(res))
     
  display_cluster(clusters)
  fig, ax = plt.subplots(figsize=(6.0, 3.0))
  try:
    graphics.plot_dendrogram(ax, tree,labels=ids,show_distance=True)
    fig.tight_layout()
  except:
    graphics.plot_dendrogram(ax, tree,show_distance=True)
    fig.tight_layout()
  

#UPGMA
file_name=sys.argv[1]
matrix,ids = ReadFile(file_name)
UPGMA(matrix,ids)