import pandas as pd
import numpy as np

trialdata = pd.read_csv('trialdata.csv')
videoinformation = pd.read_csv('videoinformation.csv')

genotypes = np.unique(trialdata['genotype'])

print(genotypes)
fly_codes = videoinformation.columns.drop(['video','trial','videodate','day','time','lightslefton','condensation','flystuck','frameshift','flicker','trackable','doubleexposed','potentialswaps','actualswaps','glitch','useable','frames','pxmm','manual_pxmm'])
fly_sex_color = []
print(fly_codes)

# for fly_code in fly_codes:
#     fly_sex_color.append(fly_code.split("."))
# print(fly_sex_color)

flynumber_to_genotype = [[]]
for index, video in videoinformation.iterrows():
    trial = video['trial']
    relevant_trialdata = trialdata[(trialdata.trial == trial)]
    for fly_code in fly_codes:
        fly_sex_color = fly_code.split(".")
        sex = fly_sex_color[0]
        color = fly_sex_color[1]
        fly_number = video[fly_code]
        row = relevant_trialdata[(relevant_trialdata.sex == sex) & (relevant_trialdata.color == color)]
        genotype = row.genotype
        flynumber_to_genotype.append([video['video'], trial, fly_number, genotype])

print(flynumber_to_genotype)

answer = pd.DataFrame(flynumber_to_genotype)

answer.to_csv("answer.csv")