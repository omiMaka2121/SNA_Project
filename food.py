import pandas as pd
import glob
import networkx as nx
import matplotlib.pyplot as plt
# Replace 'folder_path' with the actual path to your folder containing the Excel files.
csv_files = glob.glob('trialdata/*.csv')

dataframes = []
for file in csv_files:
    df = pd.read_csv(file, header=None)
    # print(df.dtypes)
    dataframes.append(df)

print(dataframes[4])


import numpy as np

# import networkx as nx
dataframe = np.array(dataframes)
print(dataframe.shape[2])
print(dataframe.shape)


graphs_edgelists = []
for d in range(56):
    row = 0
    edgelist = []
    for i in dataframe[d]:
        column = 0
        for j in i:
            if row != column:
                edgelist.append([row, column, j])
            column += 1
        row += 1
    edges_df = pd.DataFrame(edgelist, columns=['Source', 'Target', 'Weight'])
    graphs_edgelists.append(edges_df)

print(graphs_edgelists[1])


Graph_list = []
for edgelist in graphs_edgelists:
    # G = nx.from_numpy_array(np.matrix(df), create_using=nx.DiGraph)
    G = nx.from_pandas_edgelist(edgelist, source='Source', target='Target', edge_attr=True, create_using=nx.DiGraph())
    Graph_list.append(G)

print(len(Graph_list))
print(Graph_list[0].nodes)

weight = nx.get_edge_attributes(Graph_list[0], 'Weight')
print(weight)



# food 1 - vid 8th
# food 2 - vid 4th
# food 3 - vid 0th
# food 4 - vid 6th
# food 5 - vid 13th

fig, ax = plt.subplots()

eig_centrality = nx.eigenvector_centrality(Graph_list[8], weight='Weight')
list(eig_centrality.items())
df0 = pd.DataFrame(list(eig_centrality.items()), columns=['Node', 'Eigenvalue Centrality'])
plt.boxplot(df0['Eigenvalue Centrality'], widths=0.6, positions = [0])


eig_centrality = nx.eigenvector_centrality(Graph_list[4], weight='Weight')
list(eig_centrality.items())
df0 = pd.DataFrame(list(eig_centrality.items()), columns=['Node', 'Eigenvalue Centrality'])
ax.boxplot(df0['Eigenvalue Centrality'], widths=0.6, positions = [1])


eig_centrality = nx.eigenvector_centrality(Graph_list[0], weight='Weight')
list(eig_centrality.items())
df0 = pd.DataFrame(list(eig_centrality.items()), columns=['Node', 'Eigenvalue Centrality'])
ax.boxplot(df0['Eigenvalue Centrality'], widths=0.6, positions = [2])


eig_centrality = nx.eigenvector_centrality(Graph_list[6], weight='Weight')
list(eig_centrality.items())
df0 = pd.DataFrame(list(eig_centrality.items()), columns=['Node', 'Eigenvalue Centrality'])
ax.boxplot(df0['Eigenvalue Centrality'], widths=0.6, positions = [3])


eig_centrality = nx.eigenvector_centrality(Graph_list[13], weight='Weight')
list(eig_centrality.items())
df0 = pd.DataFrame(list(eig_centrality.items()), columns=['Node', 'Eigenvalue Centrality'])
ax.boxplot(df0['Eigenvalue Centrality'], widths=0.6, positions = [4])


plt.show()
