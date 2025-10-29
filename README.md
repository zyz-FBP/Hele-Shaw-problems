# Geometric local parameterization for solving Hele-Shaw-problems with surface tension
This repository contains the implementation of the methods described in our paper (https://arxiv.org/abs/2510.14088). The project focuses on solving Hele-Shaw problems with surface tension on manifolds which are identified by point clouds. 

The provided codes display the boundary evolution dynamics and reproduce the main numerical examples from the paper:

(1) FBP_circle.m - Circular case shown in Figure 3 (a) and (b) (Note that Figure 3 (b) in our paper plots the error of radius for two periods)

(2) FBP_flower.m - Perturbed circular cases shown in Figure 4 (e) and (f) (Note that Figure 4 (a)-(d) can obtained by setting different values for the parameters, D_1 and D_2), 

(3) FBP_heart.m and FBP_humanoid - Smooth closed curves shown in Figure 5.

If you use this code in your work, please cite:

@article{zhang2025geometric,
  title={Geometric local parameterization for solving Hele-Shaw problems with surface tension}, 
  author={Zhang, Zengyan and Hao, Wenrui and Harlim, John},  
  journal={arXiv preprint arXiv:2510.14088},
  year={2025}
}


