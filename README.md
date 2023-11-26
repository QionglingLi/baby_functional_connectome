# baby_functional_connectome
This repository provides data and relevant codes & toolbox used in our paper 'Development of Segregation and Integration of Functional Connectomes during the First 1000 Days'.
# overview
Content includes standalone software, source code, and demo data. Due to the large size of the analyzed data, we only provide a small portion of the data needed for validating the code. The project is structured into four parts corresponding to the major analyses in the article, including functional connectivity strength analysis, clustering coefficient & efficiency analysis, module analysis and gene expression analysis. Due to size limitation, the script relavent files and data can be found in 
# installation
Please use the “add path” method in MATLAB, "pip install" method in Python and "install.packages()" method in R to add toolboxes and scripts in the code folder. These procedures are not time-consuming.
# workflow
1. smri preprocessing
   sMRI were preprocessed with two images average, reorientation to fsl space, bias correction using N3, extract brain using skullStripping toolkit by Shifeng 2012 (https://www.nitrc.org/projects/skulltoolkit), and registration to UNC_4D_template (https://www.nitrc.org/projects/uncbcp_4d_atlas/).
2. fmri preprocessing
   
4. FCS analysis
5. graph theory measures computation
6. graph theory measures analysis
7. growth trajecotry fitting
8. gene analysis 
