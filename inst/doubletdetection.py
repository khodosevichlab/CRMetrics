def doubletdetection_py(cm, ncores = 1):
        import doubletdetection
        import io
        from contextlib import redirect_stdout
        
        clf = doubletdetection.BoostClassifier(n_jobs = ncores)
        f = io.StringIO()
        with redirect_stdout(f):
                labels = clf.fit(cm).predict()
                scores = clf.doublet_score()
        out = f.getvalue()
        return(labels, scores, out)

