def scrublet_py(cm, total_counts = None, sim_doublet_ratio = 2.0, n_neighbors = None, expected_doublet_rate = 0.1, stdev_doublet_rate = 0.02, random_state = 0, synthetic_doublet_umi_subsampling = 1.0, use_approx_neighbors = True, distance_metric = 'euclidean', get_doublet_neighbor_parents = False, min_counts = 3, min_cells = 3, min_gene_variability_pctl = 85, log_transform = False, mean_center = True, normalize_variance = True, n_prin_comps = 30, svd_solver = 'arpack'):
  import scrublet
  import io
  from contextlib import redirect_stdout
  
  tmp = scrublet.Scrublet(cm, total_counts = total_counts, sim_doublet_ratio = sim_doublet_ratio, n_neighbors = n_neighbors, expected_doublet_rate = expected_doublet_rate, stdev_doublet_rate = stdev_doublet_rate, random_state = random_state)
  f = io.StringIO()
  with redirect_stdout(f):
    doublet_scores, predicted_doublets = tmp.scrub_doublets(synthetic_doublet_umi_subsampling = synthetic_doublet_umi_subsampling, use_approx_neighbors = use_approx_neighbors, distance_metric = distance_metric, get_doublet_neighbor_parents = get_doublet_neighbor_parents, min_counts = min_counts, min_cells = min_cells, min_gene_variability_pctl = min_gene_variability_pctl, log_transform = log_transform, mean_center = mean_center, normalize_variance = normalize_variance, n_prin_comps = n_prin_comps, svd_solver = svd_solver)
  out = f.getvalue()
  return(predicted_doublets, doublet_scores, out)
