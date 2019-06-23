# Network DECODE (NetDECODE)

This repository contains code to construct sparsity-inferred co-activity networks. The main function for constructing these networks is **NetDECODE**. These networks are then pruned to construct an adaptive nearest neighbor graph using function **pruneNetwork**. Finally, activity score of each gene can be imputed as the product of pruned network by the binarized expression matrix. For more information, see *Rpackage/demo/NetDECODE_demo.Rmd* and *matlab/NetDECODE_demo.Rmd* for R and MATLAB demo files, respectively. Additionally, all C++ codes are available under *src/NetDECODE.cc*

