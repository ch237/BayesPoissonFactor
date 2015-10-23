# BayesPoissonFactor
Bayesian Poisson Tensor Factorization
This is a matlab implementation of scalable bayesian tensor factorization model, which contains both online inference
and online&parallel inference.

The code is related to our following two publications:

(1) C Hu, P Rai, C Chen, M Harding, L Carin. Scalable Bayesian Non-Negative Tensor Factorization for Massive Count Data, 
ECML-PKDD 2015, Porto, Portugal.
(2) C Hu, P Rai, L Carin. Zero-Truncated Poisson Tensor Factorization for Massive Binary Tensors,
UAI 2015, Amsterdam, The Netherlands.

Run demo_parallel.m to run the code, the data we use here is GDELT political science data.
You need to prepare the data in the form of two variables: id and xi; and store it in xi_id.mat;

id: should be a 1X4 cell, with k-th cell being a vector storing the k-th mode indices of non-zeros in the tensor  
xi: is a vector storing all the non-zeros of the tensor

You can run online inference if you set paralellFlag=0. 
If it is too slow, you can run the parallel version by setting paralellFlag=1.
If you run the parallel inference, you also need to set "M" and "cutmode". 
M is the number of cores you want to use, and cutmode is the mode along which you would like to split the data.
You can set the value of cutmode to the mode with largest dimenstionality.
