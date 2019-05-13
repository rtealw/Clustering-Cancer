'''
Spring 2019
CS321 Final Project
R. Teal Witter
Hierarchical Clustering of 
Breast Cancer Tumors and Selected Genes
'''

import sys
import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster import hierarchy
import plotly
import plotly.graph_objs as go
import plotly.io as pio

# Available PAM50 genes selected by the paper
PAM50 = np.array(["ACTR3B", "BCL2", "BLVRA", "BAG1", "ESR1", "FOXA1", "MLPH", "NAT1", "PGR", "SLC39A6", "ERBB2", "GRB7", "MMP11", "KRT14", "KRT17", "KRT5", "CDH3", "EGFR", "SFRP1", "PHGDH", "BIRC5", "MELK", "ANLN", "CCNB1", "CENPF", "CDC20", "KIF2C", "CEP55", "MKI67", "NDC80", "NUF2", "UBE2T", "TYMS", "RRM2", "UBE2C", "ORC6", "MYC"])

# PAM50 subtype for the PAM50 genes
subtypes = np.array(["Basal", "LumA", "Basal", "Basal", "Basal", "LumA", "Basal", "LumA", "LumA", "LumA", "LumA", "LumA", "LumA", "LumA", "Basal", "Basal", "Basal", "Basal", "HER2", "HER2", "HER2", "LumB", "LumB", "LumB", "Normal", "Normal", "Normal", "HER2", "HER2", "HER2", "LumB", "LumB", "LumB", "Normal", "Normal", "Normal", "HER2", "HER2", "HER2", "LumB", "LumB", "LumB", "Normal", "Normal", "Normal"])

# Order of the "flat" clusters
ORDERS = [[17,16,18,31,27,24,23,30,29,26,28,34,21,20,35,33,25,22,11,10,36,6,5,4,9,2,3,0,32,19,12,15,14,13,1,8,7],
[7,1,6,5,4,8,11,10,17,16,18,19,15,14,13,36,35,12,9,2,3,0,28,21,31,27,24,23,30,29,26,34,33,20,25,22,32],
[6,5,4,9,2,3,0,1,11,10,12,17,16,18,19,36,25,22,32,34,33,24,23,30,29,26,28,31,27,21,20,35,15,14,13,8,7],
[16,4,12,3,22,6,36,2,27,8,33,25,26,43,20,24,1,7,14,18,0,30,11,10,9,13,44,41,32,23,40,35,37,21,28,39,31,5,34,42,15,17,38,19,29],
[16,4,3,6,2,0,14,12,10,1,7,11,9,13,22,36,38,19,29,15,17,42,44,43,20,24,33,25,26,18,27,8,30,23,41,32,39,35,21,28,31,37,40,5,34],
[16,3,4,6,2,15,17,0,14,22,36,38,19,29,42,44,43,20,24,33,25,26,18,27,23,41,35,37,21,28,5,34,31,40,32,39,8,30,12,10,11,9,13,1,7]]

# The part of the CSV that contains
# protein abundance by gene, tumor
FILE_START = 34
FILE_END = 79

def read(filename):
  ''' Read from CSV file '''
  with open(filename, 'r') as csv_file:
    file = csv.reader(csv_file, delimiter=",")
    # IDs of all 45 tumors
    tumor_ids = np.array(next(file)[FILE_START:FILE_END])
    g_dict = {}
    # Read genes in PAM50
    for row in file:
      gene = row[0]
      if gene in PAM50:
        g_dict[gene] = row[FILE_START:FILE_END]
    # Check that all desired genes are in dataset 
    if sorted(g_dict.keys()) != sorted(PAM50):
      print("Some genes missing")
    # Put genes in order of PAM50
    data = []
    for gene in PAM50:
      data.append(np.array([np.log2(float(i)) for i in g_dict[gene]]))
    return np.array(data), tumor_ids

def dist_matrix(data):
  ''' Create distance matrix from data points '''
  num_pts = len(data)
  dist = np.zeros((num_pts, num_pts))
  # Compute all pairwise Euclideand distances
  for i in range(num_pts):
    for j in range(i+1, num_pts):
      dist_btwn = np.linalg.norm(data[i] - data[j])
      dist[i][j] = dist[j][i] = dist_btwn
  return num_pts, dist

def single(c1, c2, dist):
  ''' Return distance between closest pair of points in c1, c2 '''
  d = np.Inf
  for i in c1:
    for j in c2:
      if dist[i, j] < d:
        d = dist[i,j]
  return d

def complete(c1, c2, dist):
  ''' Return distance between furthest pair of points in c1, c2 '''
  d = 0
  for i in c1:
    for j in c2:
      if dist[i, j] > d:
        d = dist[i,j]
  return d

def average(c1, c2, dist):
  ''' Compute the average distance between points in c1, c2 '''
  d = 0
  for i in c1:
    for j in c2:
      d += dist[i,j] 
  return d / (len(c1) * len(c2))

def distance(c1, c2, dist, method):
  ''' Return distance according to method between c1, c2 '''
  if method.lower() == "single":
    return single(c1, c2, dist)
  if method.lower() == "complete":
    return complete(c1, c2, dist)
  if method.lower() == "average":
    return average(c1, c2, dist)
  print('Method ' + method + ' not recognized.')

def find_closest(clusters, dist, method):
  ''' Find pair of closest clusters '''
  k = len(clusters)
  clusters_dist = np.ones((k,k)) * np.Inf
  # Create matrix of distance between clusters
  for i in range(k):
    for j in range(k):
      if i != j:
        clusters_dist[i,j] = distance(clusters[i], clusters[j], dist, method)
  # Find index of minimum value in clusters_dist
  raveled = np.argmin(clusters_dist)
  i, j = np.unravel_index(raveled, clusters_dist.shape)
  return clusters[i], clusters[j]

def clusterify(n, dist, method, labels):
  ''' Execute hierarchical clustering algorithm '''
  clusters = [[cluster] for cluster in range(n)]
  # Create linkage for dendrogram
  linkage = np.empty((n-1, 4), dtype=np.float) 
  # Store which cluster we're on
  ids = [[i] for i in range(n)]
  labelled = []
  for row in range(n-1):
    ci, cj = find_closest(clusters, dist, method)   
    i, j = ids.index(ci), ids.index(cj)
    clusters.remove(ci)
    clusters.remove(cj)
    clusters.append(ci + cj)
    ids.append(ci + cj)
    linkage[row] = [i, j, distance(ci, cj, dist, method), len(ci) + len(cj)]
  return linkage

def clusterGenes(data, method):
  ''' Cluster genes by method '''
  n, dist = dist_matrix(data)
  methods = ['Single', 'Complete', 'Average']
  t= [4.2, 8, 7] # Color thresholds for dendrogram
  i = methods.index(method) # Identify the current method
  hierarchy.dendrogram(clusterify(n, dist, method, PAM50),
                       # set to t[i] for colors
                       color_threshold=0, 
                       above_threshold_color='black',
                       orientation='left',
                       labels=PAM50)
  plt.title('Hierarchical Gene Clustering by ' + method + ' Method')
  plt.ylabel('Gene')
  plt.xlabel('Distance')
  plt.savefig(method+'_gene.pdf')

def clusterTumors(data, labels, method): 
  ''' Cluster tumors by method '''
  n, dist = dist_matrix(np.transpose(data)) 
  methods = ['Single', 'Complete', 'Average']
  t = [4.2, 6, 5]
  i = methods.index(method)
  hierarchy.dendrogram(clusterify(n, dist, method, labels), 
                       color_threshold=0,
                       above_threshold_color='black',
                       leaf_rotation=270,
                       labels=labels)
  plt.title('Hierarchical Tumor Clustering by ' + method + ' Method')
  plt.ylabel('Distance')
  plt.xlabel('Tumor ID')
  plt.savefig(method+'_tumor.pdf')

def rearrange(data, method):
  ''' Rearrange data according to flat cluster orders '''
  methods = ['Single', 'Complete', 'Average'] 
  i = methods.index(method)
  row_order = np.array(ORDERS[i])   # First three correspond to gene
  col_order = np.array(ORDERS[i+3]) # Last three correspond to tumor
  new_data = data[:,col_order]
  new_data = new_data[row_order,:]
  return new_data, col_order, row_order

def printPAM(col_order):
  ''' Print shorter PAM50 subtypes for easier viewing '''
  old = ['Basal', 'LumB', 'HER2', 'Normal', 'LumA']
  new = ['Ba', 'LB', 'HE', 'No', 'LA']
  for pam in subtypes[col_order]:
    print(new[old.index(pam)], end=" ")
  print()

def cluster(data, tumor_ids, method):
  ''' Cluster data according to method '''
  clusterGenes(data, method)
  clusterTumors(data, tumor_ids, method)
  # Rearrange data so heatmap corresponds to clustering
  new_data, col_order, row_order = rearrange(data, method)
  printPAM(col_order)
  trace = go.Heatmap(z=new_data,
                     x=tumor_ids[col_order],
                     y=PAM50[row_order[::-1]],
                     colorscale=[[0.0,'rgb(12,0,255)'],
                                 [0.5, 'rgb(255,255,255)'],
                                 [1.0, 'rgb(250,13,0)']])
  layout = go.Layout(
    title=method+' Heatmap',
    xaxis = dict(ticks='', nticks=45),
    yaxis = dict(ticks='', nticks=37 ),
    height= 700
  )
  fig = go.Figure([trace], layout=layout)
  pio.write_image(fig, method+'_heatmap.pdf')
  #plotly.offline.plot([trace], filename=method+'_heatmap.html')
  

if len(sys.argv) != 2:        
  print("use: python3 cluster.py <filename.csv>")
else:
  data, tumor_ids = read(sys.argv[1])
  for method in ['Single', 'Complete', 'Average']:
    cluster(data, tumor_ids, method) 




















