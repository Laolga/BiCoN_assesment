import time
import pandas as pd
from bicon import *

gene_expression_files = ["breast_expr_all3.csv", "breast_expr_ab.csv", "breast_expr_ab_ni.csv", "breast_expr_basal.csv"]

for expr_name in gene_expression_files:
    print("{0} is processing".format(expr_name))
    prefix = '/nfs/home/users/olgala/scratch/data_processing/tcga_processed/'
    path_expr = prefix + expr_name
    path_net = '/nfs/home/users/olgala/scratch/biogrid/biogrid_net.tsv'
    GE, G, labels, rev = data_preprocessing(path_expr, path_net, size = 3000, log2=True)
    results = []
    print("{0} analysis starts".format(expr_name))
    for sim in range(10):
        st = time.time()
        model = BiCoN(GE, G, L_g_min=10, L_g_max=25)
        solution, scores = model.run_search(ls=True, n_proc = 10)
        end = time.time()
        patients1 = [str(labels[x]) for x in solution[1][0]]
        patients1 = '|'.join(patients1)

        patients2 = [str(labels[x]) for x in solution[1][1]]
        patients2 = "|".join(patients2)

        genes1 = [str(labels[x]) for x in solution[0][0]]
        genes1 = "|".join(genes1)
        sc = scores[1][-1]

        genes2 = [str(labels[x]) for x in solution[0][1]]
        genes2 = "|".join(genes2)

        cols = ["score", "time", "patients1", "patients2",
                "genes1", "genes2"]
        sol = [sc, round(end - st), patients1, patients2, genes1, genes2]

        results.append(sol)
        results_df = pd.DataFrame(results, columns=cols)
        results_df.to_csv("../results/" + expr_name, index=False)

