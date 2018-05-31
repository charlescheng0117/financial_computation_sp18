Homework 3
===

Link: http://homepage.ntu.edu.tw/~jryanwang/course/Financial%20Computation%20or%20Financial%20Engineering%20(graduate%20level)/Homework_3.pdf

Basic Requirement
---
- Rainbow Option: Max( Max(S1, S2, ..., Sn) - K, 0 )
- Cholesky Decomposition Method: Algorithm is on FE_Ch05, Page 5-4


Bonus 1
---
- Antithetic variance approach: satisfying the zero mean and the symmetric feature of the standard normal distribution. FE_Ch05, Page 5-8, (i)
- Moment matching: matching the first two memonts of the standard normal distribution. FE_Ch05, Page 5-8, (ii)
    - Draw random samples z1, z2, ..., zn first. Suppose mean = m, s.d. = s.
    - Define: 
        - yi = (zi - m) / s => get y1, y2, ..., yn
        - mean = 0, s.d. = 1 now

Bonus 2
---
- Wang (2008), “Variance Reduction for Multivariate Monte Carlo Simulation,” Journal of Derivatives 16, pp. 7–28.
- The New Method
    - Step 1: Generate independent standard normal distributed random samples for each underlying asset and obtain a matrix of random samples.
    - Step 2: Calculate the variance-covariance matrix C between zj^~ = zj - uj, zk^~ = zk - uk, where uk, uj are the sample means of zj, jk.
    - Step 3: Based on the covariance matrix C~, perform the Cholesky decomposition C~ = ÃTÃ to obtain the corresponding linear transformation Ã. Find the inverse matrix
Ã^-1
    - Step 4: : Since the variance-covariance matrix of [z˜1
z˜
2 … z˜
N] is not exactly the identity matrix, it can be viewed
as a group of correlated normally distributed random samples
obtained from [z1
´z2
´ … zN
´] × Ã, where [z1
´z2
´ … zN
´]
are truly independent, standard normally distributed random
samples. By applying the inverse Cholesky decomposition
matrix Ã–1, the matrix of truly independent standard normal
distributed random samples Z´ can be obtained by
Z´= [z1
´z2
´ … zN
´] = [z˜ 1z˜ 2 … z˜N] × Ã–1 =Z˜ × Ã–1 
