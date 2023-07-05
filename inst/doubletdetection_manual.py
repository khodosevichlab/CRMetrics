print("Importing libraries")

# Setup
import doubletdetection
from scipy import io
import pandas as pd
import glob

flist = glob.glob("data.path*.mtx")

for x in flist:
    print("Creating "+x+".method")
    
    # Load
    cm = io.mmread(x)

    # Run
    clf = doubletdetection.BoostClassifier(boost_rate = boost_rateX, clustering_algorithm = clustering_algorithmX, clustering_kwargs = clustering_kwargsX, n_components = n_componentsX, n_iters = n_itersX, n_jobs = n_jobsX, n_top_var_genes = n_top_var_genesX, normalizer = normalizerX, pseudocount = pseudocountX, random_state = random_stateX, replace = replaceX, standard_scaling = standard_scalingX)
    labels = clf.fit(cm).predict(p_thresh = p_threshX, voter_thresh = voter_threshX)
    scores = clf.doublet_score()

    # Save
    pd.DataFrame(labels, scores).to_csv(x+".method")
