
import time
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib
matplotlib.use('agg')
flatten = lambda l: [item for sublist in l for item in sublist]
import seaborn as sns; sns.set(color_codes=True)
import bicon


from bicon import *

path_expr, path_net = '../data_processing/GSE30219/lung_expr_nonorm.csv', '../biogrid/biogrid_net.tsv'
GE, G, labels, rev = data_preprocessing(path_expr, path_net, size = 3000, log2=True)
n,m = GE.shape
new_pat = np.arange(n,n+m)
new_genes = np.arange(n)


results = []

props = np.arange(0,1, step = 0.1)

for i in props:
    print(i)
    new_expr = GE.copy().values
    size = int((i)*m)
    for l in range(n):
        col_id = list(np.random.choice(np.arange(m), size=size, replace=False))
        v = new_expr[l, col_id]
        np.random.shuffle(v)
        new_expr[l, col_id] = v
    new_expr = pd.DataFrame(new_expr,index = new_genes , columns= new_pat)
    for runs in range(10):
        st = time.time()
        model = BiCoN(GE, G, L_g_min = 10, L_g_max = 25)
        solution, scores = model.run_search(ls=True)
        end = time.time()
        patients1 = [str(labels[x]) for x in solution[1][0]]
        patients1 = "|".join(patients1)

        patients2 = [str(labels[x]) for x in solution[1][1]]
        patients2 = "|".join(patients2)

        genes1 = [str(labels[x]) for x in solution[0][0]]
        genes1 = "|".join(genes1)

        genes2 = [str(labels[x]) for x in solution[0][1]]
        genes2 = "|".join(genes2)
        sc = scores[1][-1]


        cols = ["prop", "algorithm", "score","time","size_min","size_max","patients1","patients2","genes1","genes2"]
        complexity = n*m*len(G.edges)
        sol = [i,"BiCoN",sc,round(end-st),10,25,patients1,patients2,genes1,genes2]

        results.append(sol)
        results_df = pd.DataFrame(results,columns = cols)
        results_df.to_csv("results/robustness.csv", index = False)
