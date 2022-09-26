def doubletdetection_py(cm, boost_rate = 0.25, clustering_algorithm = "phenograph", clustering_kwargs = None, n_components = 30, n_iters = 10, n_jobs = 1, n_top_var_genes = 10000, normalizer = None, pseudocount = 0.1, random_state = 0, replace = False, standard_scaling = False, p_thresh = 1e-7, voter_tresh = 0.9):
        import doubletdetection
        import io
        from contextlib import redirect_stdout
        
        clf = doubletdetection.BoostClassifier(boost_rate = boost_rate, clustering_algorithm = clustering_algorithm, clustering_kwargs = clustering_kwargs, n_components = n_components, n_iters = n_iters, n_jobs = n_jobs, n_top_var_genes = n_top_var_genes, normalizer = normalizer, pseudocount = pseudocount, random_state = random_state, replace = replace, standard_scaling = standard_scaling)
        f = io.StringIO()
        with redirect_stdout(f):
                labels = clf.fit(cm).predict(p_thresh = p_thresh, voter_tresh = voter_tresh)
                scores = clf.doublet_score()
        out = f.getvalue()
        return(labels, scores, out)

