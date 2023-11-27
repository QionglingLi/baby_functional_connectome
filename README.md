# baby_functional_connectome
This repository provides data and relevant codes & toolbox used in our paper 'Development of Segregation and Integration of Functional Connectomes during the First 1000 Days'.
# overview
Content includes standalone software, source code, and demo data. Due to the large size of the analyzed data, we only provide a small portion of the data needed for validating the code. The project is structured into four parts corresponding to the major analyses in the article, including functional connectivity strength analysis, clustering coefficient & efficiency analysis, module analysis and gene expression analysis. Due to size limitation, the script relavent files and data can be found in https://pan.bnu.edu.cn/l/N1uFm7
# installation
Please use the “add path” method in MATLAB, "pip install" method in Python and "install.packages()" method in R to add toolboxes and scripts in the code folder. These procedures are not time-consuming.
# workflow
1. **smri preprocessing**
   sMRI were preprocessed with following steps: average two images, reorient to fsl space, bias correction using N3, extract brain using [skullStripping toolkit](https://www.nitrc.org/projects/skulltoolkit) by Shifeng, and registration to [UNC_4D_Volume_template](https://www.nitrc.org/projects/uncbcp_4d_atlas/). Details can be found in the manuscript and the [script](https://github.com/QionglingLi/baby_functional_connectome/blob/main/code/sMRI_Preprocessing.sh).
2. **fmri preprocessing**
   fMRI were preprocessed with following steps: drop first ten TRs and reorient, motion correction with SBRef, distortion correction, SBRef to anatomical registration, registration to BCP age-specific template and then to common template, smoothing, detrend & nuisance regressing, bandpass filtering, and scrubbing. Details can be found in the manuscript and the [script](https://github.com/QionglingLi/baby_functional_connectome/blob/main/code/fMRI_Preprocessing.sh).    
3. **FCS analysis**
   FCS were computed at the voxel level. The gray matter mask were predefined by applying a threshold to the [gray matter probability template of 6 months](https://github.com/QionglingLi/baby_functional_connectome/blob/main/data/BCP-06M-GM_2mm.nii.gz). The FCS (also within each distance bin) computation for each individual and the hub identification were computed using codes [here](https://github.com/QionglingLi/baby_functional_connectome/blob/main/code/1_FCS%20analysis.m). In order to explore the distribution of hubs across different functional systems, we mapped the Yeo's atlas on the 6-month infant brain. We first nonlinearly registered the MNI152 template to the BCP 24-month template using [ANTs SyN algorithm](https://github.com/ANTsX/ANTs). And then, we linearly registered it to the 6-month template. The combined deformation was applied to Yeo's 7 network atlas in MNI152 space. Yeo's 7 network on 6-month infant brain can be found [here](https://github.com/QionglingLi/baby_functional_connectome/blob/main/data/Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask_BCP24_to_BCP06_2mm.nii.gz).
4. **graph theory measures computation**
   The graph theoretical measures on voxel-wise brain networks were executed using the [Parallel Graph-theoretical Analysis (PAGANI) toolbox](https://www.nitrc.org/projects/pagani_toolkit/).
5. **graph theory measures analysis**
    The principal development axis of graph theoretical measures were analyzed using codes [here](https://github.com/QionglingLi/baby_functional_connectome/blob/main/code/2_GraphTheory_Analysis.m).
6. **growth trajectory fitting**
    The developmental trajectories were fitted by generalized additive mixed model (GAMM) in R. The dependent package is mgcv (version 1.8-42). The Gaussian random field (GRF) correction for the voxel level analysis were using the toolbox [SeeCAT](https://www.nitrc.org/projects/seecat/). 
7. **gene analysis**
    Gene analysis and figure codes are all from Zhilei Xu's work "Meta-connectomic analysis maps consistent, reproducible, and transcriptionally relevant functional connectome hubs in the human brain" pusblised in the Communications Biology.
The code link is [here](https://github.com/zhileixu/FunctionalConnectomeHubs/tree/main/Figure6/Figure6a) & [here](https://github.com/zhileixu/FunctionalConnectomeHubs/tree/main/Figure6/Figure6b).
8. **visualization**
    The cortical surface visualization in this work were using the [BrainNet Viewer](http://www.nitrc.org/projects/bnv/) software. The cortical surface file is [here](https://github.com/QionglingLi/baby_functional_connectome/blob/main/data/BCP_06month.nv).
   
# citation
please cite us as follows:
