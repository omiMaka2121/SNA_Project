import pandas as pd
import glob
import networkx as nx

# Replace 'folder_path' with the actual path to your folder containing the Excel files.
csv_files = glob.glob('trialdata/*.csv')

dataframes = []
for file in csv_files:
    df = pd.read_csv(file, header=None)
    # print(df.dtypes)
    dataframes.append(df)


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

# weight = nx.get_edge_attributes(Graph_list[0], 'Weight')
# print(weight)

###

answers = pd.read_csv('answer.csv')
answers['fly_number'] = answers['fly_number'] - 1
print(answers)

import os
# specify the path to the folder
folder_path = 'trialdata'
# get a list of all filenames in the folder
filenames = os.listdir(folder_path)
print(filenames)  # prints a list of filenames in the folder
filenames = [s.replace('.csv', '') for s in filenames]
print(filenames)

filtered_trials=[]
for i in filenames:
    file_current=i
    filtered_df = answers.loc[answers['video'] == file_current]
    filtered_trials.append(filtered_df)
    print(filtered_trials[0])
# Print the dictionary
for i in filtered_trials:
    pass

#Final Renaming of nodes
Renamed_Graphs=[]
import pandas as pd
j=0
for i in filenames:
    file_current=i
    filtered_df = answers.loc[answers['video'] == file_current]
    filtered_df=pd.DataFrame(filtered_df)
    print(filtered_df)
    renaming_dict = {k: v for k, v in zip(filtered_df['fly_number'], filtered_df['genotype'])}
    print(renaming_dict)
    H= nx.relabel_nodes(Graph_list[j], renaming_dict)
    Renamed_Graphs.append(H)
    print(j)
    j=j+1
    if(j==56):
        break

#Getting the genotype list
Genotype_list=answers['genotype'].tolist()
#print(Genotype_list)
Genotype_list_unique = list(set(Genotype_list))
# print the new list without duplicates
print(Genotype_list_unique)

Genotypes_Eigenvalues=[]
j=0
for i in Renamed_Graphs:
    eig_centrality = nx.eigenvector_centrality(i, weight='Weight')
    # create a dataframe from the eigenvalue centrality dictionary
    if(j==0):
        Centrality_df = pd.DataFrame(list(eig_centrality.items()), columns=['Node', 'Eigenvalue Centrality'])
    else :
        new_row = pd.DataFrame(list(eig_centrality.items()), columns=['Node', 'Eigenvalue Centrality'])
        Centrality_df=pd.concat([Centrality_df, new_row])
    j=j+1
print(Centrality_df.head(10))

filtered_Cent_df = Centrality_df[Centrality_df['Node'] == '732x40' ]
print(filtered_Cent_df)

import matplotlib.pyplot as plt
# plt.boxplot(filtered_Cent_df['Eigenvalue Centrality'])

# add labels and title
plt.xlabel('Group')
plt.ylabel('Value')
plt.title('Box Plot Example')

# show the plot
# plt.show()

import matplotlib.pyplot as plt
import numpy as np

# create some sample data
np.random.seed(1234)
data = [np.random.normal(0, std, 100) for std in range(1, 4)]
print(data)
# create a box plot
plt.boxplot(data)
# add labels and title
plt.xlabel('Group')
plt.ylabel('Value')
plt.title('Box Plot Example')
# show the plot
# plt.show()

fig, ax = plt.subplots()

for i in range(len(Genotype_list_unique)):
    filtered_Cent_df = Centrality_df[Centrality_df['Node'] == Genotype_list_unique[i] ]
    ax.boxplot(filtered_Cent_df['Eigenvalue Centrality'], positions=[i], widths=0.6)
    # ax.set_xticks(Genotype_list_unique[i])
    # add labels and title
    plt.xlabel(Genotype_list_unique[i])
    plt.ylabel('Value')
    plt.title('Box Plot Example')

    # show the plot

    print(i)

# plt.show()

# import matplotlib.pyplot as plt
import numpy as np

# Generate some random data for the box plots

# Add labels and title to the plot
ax.set_xticks([i for i in range(20)])
ax.set_xticklabels(Genotype_list_unique)
ax.set_ylabel('Values')
ax.set_title('Box plots of multiple data sets')

# Show the plot
plt.xticks(rotation=30)
plt.show()