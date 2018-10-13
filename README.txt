Thank you for downloading this MATLAB package for Joint Matrix Factorization. This package contains the following files.

This folder contains Matlab code related to the following paper:

Lihua Zhang, Shihua Zhang. A Unified Joint Matrix Factorization Framework for Data Integration. Preprint arXiv:1707.08183. 
(https://arxiv.org/abs/1707.08183)

The folder was last updated by Lihua Zhang on July 31, 2017.

Please look `main_demo.m' at first and try to execute examples to learn how to use this software. Please send bug reports, comments, or questions to zsh@amss.ac.cn.

=================================
DESCRIPTION
=================================

The codes including three part: simulation part, main function part, prediction part
%%%%%%% simulation part
JMF_synthetic_dataset1.m
JMF_synthetic_dataset2.m
JMF_synthetic_dataset3.m
JMF_synthetic_dataset4.m
These four codes used for generating synthetic datasets.

%%%%%%% mian function part
The main function is named JMF.m, which includes many parameters used to select update rules and stop criterions. More details are illustrated in the function.
MUR_updateW.m, MUR_updateH.m
PG_updateW.m, PG_updateH.m
Ne_updateW.m, Ne_updateH.m
PANLS_updateW.m, PANLS_updateH.m
There codes used for updating W and H, which are invoked by JMF.m, JMF_prediction_L.m, and JMF_prediction_R.m.

%%%%%%% prediction part
Two prediction codes are included in this part, which are named by JMF_prediction_L.m and JMF_prediction_R.m.

%%%%%%% Other codes (compute_diff.m, compute_sumHRt.m, compute_sumtheta_record.m, compute_XHt_HHt.m, get_gradW.m, get_gradH.m, StopCriterion_rule1.m, and StopCriterion_rule2.m) are invoked by many codes mentioned above.






