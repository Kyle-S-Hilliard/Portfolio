'''
BASIC SINGLE CELL ANALYSIS SCRIPT
by Josh Wu
4 June, 2019

Relies heavily on the Scanpy Python module developed by the Theis Lab
Read more about Scanpy at https://scanpy.readthedocs.io/en/latest/index.html

Contains analysis of Mouse intestine samples for Kate Walton

In progress ---
Moving to encapsulate parameters and relevant functions using class sca_set()
'''
import sys
sys.path.insert(0,'C:/Users/Josh/Desktop/sca/tools')

from sca_run import *
#from tools.pipelines import *

figdir = './figures_012920/'
an_run = sca_run()
#############################################################################
## Change this to point toward your mount location for our MiStorage share ##
#############################################################################
an_run.storage_mount_point = 'Z:/'

## IDs of samples as represented in the metadata table
an_run.sample_list = ['186-1','186-2','186-3']

# ## List of interesting genes
an_run.add_gene_list(markers = ['Vim','Cdh2','Cdh1','Vil1','Shh','Ihh','Acta2','Tagln','Acta1','Des','Smtn',
								'Myh11','Foxl1','Foxf1','Pdgfra','Cd34','Gli1','Pdpn','Pecam1','Kdr','Cdh5',
								'Lyve1','Flt4','Cspg4','Pdgfrb','Mcam','Nes','Ptprc','Adgre1','Sirpa','Mertk',
								'Cd68','Itgam','Itgax','Itgae','Cd19','Il7r','Ms4a1','Cd93','Tcf3','Pax5',
								'Ebf1','Pou2f2','Cd4','Cd8a','Tubb3','Prph','Ngfr','Ret','S100b','Gfap','Wt1',
								'Msln','Kit','Calb2','Ano1']#,'HLA-B','Hla-b'],#'Mhcii'
					 label='general_list')

an_run.add_gene_list(markers = ['Shh','Ptch1','Ihh','Gli1','Gli2','Gli3','Smo','Spop','Hhip','Sufu'], #ptc1
					 label='HH_signaling')

an_run.add_gene_list(markers = ['Bmp2','Bmp4','Bmp5','Bmp7','Nog','Chrd','Grem1','Chrdl1'],
					 label='bmp_list')

an_run.add_gene_list(markers = ['Dll1','Dll4','Notch1','Notch2','Notch3','Notch4','Srrt','Jag1','Jag2','Mmp9',
								'Numb','Mfng','Lfng','Rfng'],
					 label='notch_list')

## Parameters used to filter the data - Mainly used to get rid of bad cells
an_run.set_filter_params(min_cells = 0, # Filter out cells 
						 min_genes = 500, # Filter out cells with fewer genes to remove dead cells
						 max_genes = 7500, # Filter out cells with more genes to remove most doublets
						 max_counts = 60000, # Filter out cells with more UMIs to catch a few remaining doublets
						 max_mito = 0.1) # Filter out cells with high mitochondrial gene content

## Parameters used for initial clustering analysis
an_run.set_analysis_params(n_neighbors = 30, # Size of the local neighborhood used for manifold approximation
						   n_pcs = 20, # Number of principle components to use in construction of neighborhood graph
						   spread = 1, # In combination with min_dist determines how clumped embedded points are
						   min_dist = 0.4, # Minimum distance between points on the umap graph
						   resolution = 0.4) # High resolution attempts to increases # of clusters identified

an_run.set_plot_params(size = 5,
					   umap_obs = ['louvain','sampleName'],
					   exp_grouping = ['louvain'],
					   umap_feature_color = 'yellow_blue')
					   #umap_categorical_color = ['#33EB33','#000766','#1A6DCC','#005200','#009F00','#ED8D93'])

## Basic pipeline for analysis - will filter data, process, cluster, etc. and output relevant figures
an_run.pipe_basic(figdir,load_save='adata_save.p')

## If you find some interesting clusters that you want to "zoom in" on and recluster, you can use the following code
# New analysis parameters for the subset of parameters
analysis_params_ext = dict(n_neighbors = 9,
						n_pcs = 10,
						spread = 1,
						min_dist = 0.4,
						resolution = 0.4)

an_run.umap_cluster_color = 'default'
an_run.size=160
#an_run.gene_dict['new_list']['groupby_positions'] = None
umap_categorical_color='default'

#an_run.pipe_ext(analysis_params_ext, figdir=figdir, extracted=['5'], label='cluster_5/', load_save='adata_save.p')