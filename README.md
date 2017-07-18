# pyNBS
This is the python implementation of Network-Based Stratification (NBS, https://www.nature.com/nmeth/journal/v10/n11/full/nmeth.2651.html) with some changes (see below). 

pyNBS.py(./pyNBS.py) contains the functions for network propagation and stratification. 

pyNBS_cookbook.ipynb shows an example to stratify TCGA Lung Adenocarcinoma patients from a Mutation Annotation Format (MAF) files.

label2coxph.R is the code for survival analysis.

1. Network propagation is used to compute the pairwise similarities among tumor genetic alteration profiles (0 = wild type; 1 = altered) within a refernce network. Each tumor genetic profile is propagated across this network based on a random walk model with a restart probability of 0.2-0.5. After convergence, the score of each gene represents its network proximity to genetic alterations.

2. The top 50-100 principal components of these scores, representing the tumor’s network-transformed profile are analyzed using spectral clustering (k-means and hierarchical clustering are also surpported). This method first constructs a similarity graph on all pairs of tumors, where each tumor is connected to the k others with shortest Euclidean distance. k should be chosen to ensure the similarity graph to be connected. Next this graph is analyzed to partition tumors into subtypes at different resolutions (number of subtypes n = [2..10]).

3. We use the “coxph” package in the R statistics software to fit Cox proportional hazard models. P-values are calculated by the log likelihood ratio test. 

4. For each subtype, we define a set of ‘signature genes’ as those that have higher network-transformed scores in that subtype than others (t-test) and, among these, are more frequently altered in that subtype (Fisher Exact Test). To identify subnetworks, this set is expanded to include ‘intermediate genes’ with relatively high network-transformed score (t-test) that lie on the shortest paths between each pair of signature genes. The union of the signature and intermediate genes is used to induce a subnetwork impacted in that subtype.
