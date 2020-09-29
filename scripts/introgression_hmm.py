#!/usr/bin/env python3

"""
Matt Gibson
Aug. 2019
Indiana University Departement of Biology
Moyle Lab

Viterbi for Galapagos introgression analysis


Arg 1: input file containing probabilities(generated with R script)
Arg 2: output annotation file name
Arg 3: 1 to write output 0 to not
"""



import sys
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt



"""
getTrace() performs the backtracing step in the viterbi algorithm to return a state path
"""
def getTrace(matrix, trace, states):
    lastCol = matrix.shape[1]
    last = list(matrix[:,lastCol-1])
    maxx = max(last)
    startState = last.index(maxx)

    path = []    
    currentState = startState
    for x in range(lastCol, 1, -1):
        col = list(trace[:,x-1])
        prevState = col[int(currentState)].decode("utf-8").split(',')[0]
        path.append(states[int(prevState)])
        currentState = prevState
    return(path[::-1])

"""
max_over_k() is used to find the max probability over all states in viterbi
"""
def max_over_k(prevColumn, trans, BB, BB1, BB2, state, states, i):
    numKs = len(prevColumn)
    possibleValues = []
    for k in range(numKs):
        T = float(prevColumn[k])
        AlookupValue = states[k]+state
        A = float(trans[AlookupValue])
        if state == '+':
            possibleValues.append(T + math.log10(A)+BB1)
        elif state == '-':
            possibleValues.append(T+math.log10(A)+BB)
        elif state == 'h':
            possibleValues.append(T+math.log10(A)+BB2)
    maxVal = max(possibleValues)
    k_ind = possibleValues.index(maxVal)
    argMax = str(k_ind) +',' + str(i-1)
    return(max(possibleValues), argMax)


def readData(file):
    return(pd.read_csv(file, sep = '\t'))




##############################################
#Read in observations (generated in R)
##############################################
print('READING FILE...')
df = readData(sys.argv[1])

p_114 = list(df['P_114'])
p_115 = list(df['P_115'])
p_het = list(df['P_het'])


##############################################
#Initial probabilities (uniform)
##############################################
pi = {
    '+':float(1/3),
    '-':float(1/3),
    'h':float(1/3)
}


##############################################
#Three possible states
##############################################
states = ['+', '-', 'h']


##############################################
#Transition probabilities
##############################################
admix_time = 10
r = 0.00000002*100000*admix_time
admix_proportion = 0.48

#Define transition probabilities as functions of admix_time, t, and admix_proportion
plus_plus = ((1-r)+r)*(admix_proportion)
plus_minus = r*admix_proportion
minus_minus = ((1-r)+r)*(1- admix_proportion)
minus_plus = r*admix_proportion

plusminus_plusminus = ((1-r)+r)*(1- admix_proportion)
plus_plusminus = r*admix_proportion
minus_plusminus = r*admix_proportion

#Working with windows, so we dont take into account map distance between symbols
trans = {'++':plus_plus, '+-':plus_minus, '-+':minus_plus, '--':minus_minus,
        'hh':plusminus_plusminus, '+h':plus_plusminus, '-h':minus_plusminus, 'h+':plus_plusminus, 'h-':minus_plusminus}
print(trans)





##############################################
#Start Viterbi
##############################################

print("INITIALIZING MATRIX...")
#Initialize DP matrix
i = len(states)
j = len(p_115)
print(i,j)
alg = 0
iL = list(range(0,i))
jL = list(range(0,j+1))
matrix = np.zeros(shape=(i,j+1), dtype=np.float64)
trace = np.zeros(shape=(i,j+1),dtype='<S16')
matrix[:,0] = iL
rows = matrix.shape[0]
cols = matrix.shape[1]
i = 0
for r in range(rows):
    currentState = states[i]
    iniLookupValue = currentState
    if i == 0:
        matrix[r,0] = 1
    else:
        BB = float(p_115[0])
        BB1 = float(p_114[0])
        BB2 = float(p_het[0])
        if currentState == '+':
            matrix[r,0] = math.log10(float(pi[iniLookupValue])) + float((BB1))
        elif currentState == '-':
            matrix[r,0] = math.log10(float(pi[iniLookupValue])) + float(BB)
        elif currentState == 'h':
            matrix[r,0] = math.log10(float(pi[iniLookupValue])) + float(BB2)
    i += 1

print("RUNNING VITERBI...")
for i in range(1, cols):
    for j in range(rows):
        prevColumn = matrix[:,i-1]
        maximum, arg = max_over_k(prevColumn, trans, float(p_115[i-1]), float(p_114[i-1]), float(p_het[i-1]), states[j], states, i)
        matrix[j,i] = maximum
        trace[j,i] = arg
print("BACKTRACING...")
path = getTrace(matrix, trace, states)
a = ''.join(path)
print(a)
df['annt'] = list(a)
df['pos'] = list(range(0, len(df['dxy'])))


###Optional plotting for quick debugging###
###########################################
"""
positive = df[df['annt'] == '+']
negative = df[df['annt'] == '-']
het = df[df['annt'] == 'h']

fig, ax = plt.subplots()

ax.scatter(positive['pos'], positive['dxy'], c = 'red')
ax.scatter(negative['pos'], negative['dxy'], c = 'blue')
ax.scatter(het['pos'], het['dxy'], c = 'purple')

plt.show()
"""
###########################################
###########################################


#Write output file if third arg is 1
if sys.argv[3] == '1':
    print(df)
    df.to_csv(sys.argv[2], sep = '\t')

