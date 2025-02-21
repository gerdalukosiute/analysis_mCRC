import celltypist
import scanpy as sc
from celltypist import models

models.download_models(model = 'Adult_Mouse_Gut.pkl')
predictions = celltypist.annotate(filename = '/Users/gerdalukosiute/Downloads/Thesis/seurat_expression.csv', majority_voting = True, model='Adult_Mouse_Gut.pkl')
predictions.to_table(folder = '/Users/gerdalukosiute/Downloads/Thesis', prefix = '')
predictions.to_plots(folder = '/Users/gerdalukosiute/Downloads/Thesis', prefix = '')