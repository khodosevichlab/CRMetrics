def scrublet_py(cm):
  import scrublet
  import io
  from contextlib import redirect_stdout
  
  tmp = scrublet.Scrublet(cm)
  f = io.StringIO()
  with redirect_stdout(f):
    doublet_scores, predicted_doublets = tmp.scrub_doublets()
  out = f.getvalue()
  return(predicted_doublets, doublet_scores, out)
