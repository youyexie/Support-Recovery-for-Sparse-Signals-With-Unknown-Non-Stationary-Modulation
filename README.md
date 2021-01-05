# Support-Recovery-for-Sparse-Signals-With-Unknown-Non-Stationary-Modulation
Code to plot the figures in the IEEE Transaction on Sigal Processing (TSP) paper "[Support Recovery for Sparse Signals With Unknown Non-Stationary Modulation](https://ieeexplore.ieee.org/abstract/document/9007495)"

# Abstract
The problem of estimating a sparse signal from low dimensional noisy observations arises in many applications, including super resolution, signal deconvolution, and radar imaging. In this paper, we consider a sparse signal model with non-stationary modulations, in which each dictionary atom contributing to the observations undergoes an unknown, distinct modulation. By applying the lifting technique, under the assumption that the modulating signals live in a common subspace, we recast this sparse recovery and non-stationary blind demodulation problem as the recovery of a column-wise sparse matrix from structured linear observations, and propose to solve it via block ℓ 1 -norm regularized quadratic minimization. Due to observation noise, the sparse signal and modulation process cannot be recovered exactly. Instead, we aim to recover the sparse support of the ground truth signal and bound the recovery errors of the signal's non-zero components and the modulation process. In particular, we derive sufficient conditions on the sample complexity and regularization parameter for exact support recovery and bound the recovery error on the support. Numerical simulations verify and support our theoretical findings, and we demonstrate the effectiveness of our model in the application of single molecule imaging.

# Tested on 
- Matlab R2017b with [CVX toolbox](http://cvxr.com/cvx/)

# Citation
If you use our method and/or codes, please cite our paper

```
@article{xie2020support,
  title={Support Recovery for Sparse Signals With Unknown Non-Stationary Modulation},
  author={Xie, Youye and Wakin, Michael B and Tang, Gongguo},
  journal={IEEE Transactions on Signal Processing},
  volume={68},
  pages={1884--1896},
  year={2020},
  publisher={IEEE}
}
```
