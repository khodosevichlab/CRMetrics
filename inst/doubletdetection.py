def doubletdetection_py(cm):
        import doubletdetection
        import io
        from contextlib import redirect_stdout
        
        clf = doubletdetection.BoostClassifier()
        f = io.StringIO()
        with redirect_stdout(f):
                labels = clf.fit(cm).predict()
                scores = clf.doublet_score()
        out = f.getvalue()
        return(labels, scores, out)

