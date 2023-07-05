print("Importing libraries")

import scrublet
from scipy import io
import pandas as pd
import glob

flist = glob.glob("data.path*.mtx")

for x in flist:
    print("Creating "+x+".method")
    
    # Load
    cm = io.mmread(x)
    
    # Run
    tmp = scrublet.Scrublet(cm, total_counts = total_countsX, sim_doublet_ratio = sim_doublet_ratioX, n_neighbors = n_neighborsX, expected_doublet_rate = expected_doublet_rateX, stdev_doublet_rate = stdev_doublet_rateX, random_state = random_stateX)
    doublet_scores, predicted_doublets = tmp.scrub_doublets(synthetic_doublet_umi_subsampling = synthetic_doublet_umi_subsamplingX, use_approx_neighbors = use_approx_neighborsX, distance_metric = distance_metricX, get_doublet_neighbor_parents = get_doublet_neighbor_parentsX, min_counts = min_countsX, min_cells = min_cellsX, min_gene_variability_pctl = min_gene_variability_pctlX, log_transform = log_transformX, mean_center = mean_centerX, normalize_variance = normalize_varianceX, n_prin_comps = n_prin_compsX, svd_solver = svd_solverX)

    # Save
    pd.DataFrame(predicted_doublets, doublet_scores).to_csv(x+".method")
