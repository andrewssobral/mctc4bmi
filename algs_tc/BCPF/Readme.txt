Matlab code of Bayesian CP Factorization for Tensor Completion
(Written by Qibin Zhao 2014)


To run the code:
1. Change Matlab work directory to "/BCPF_Toolbox_QZhao/".
2. Run  "loadpah" code to add the current folder and subfolders into Matlab path searching list.
3. Open and run the demo files. 


We provide two demo codes:
I.  DemoBayesCP.m:        Demonstration on synthesic data
II. DemoBayesCP_Image.m   Demonstration for image completion


The package includes four algorithms:
1. BCPF.m           BCPF for fully observed tensor
2. BCPF_TC.m        BCPF for incomplete tensor 
3. BCPF_IC.m        BCPF for image completion
4. BCPF_MP.m        BCPF using mixture priors for image completion


In this package, we used the tensor toolbox 2.5, which is downloaded from (http://www.sandia.gov/~tgkolda/TensorToolbox)

The tools for visualization of tensor with voxels is from Tensorlab (http://www.tensorlab.net/)