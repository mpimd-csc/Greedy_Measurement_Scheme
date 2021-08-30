## Adaptive selection of interpolation points for data driven modelling

This repository contains the MATLAB implementation for the adaptive construction of a data driven realization from frequency data or time domain data. The algorithm consists in the adaptive selection of interpolation points used in the construction of data driven models using the Loewner framework [1]. In each step, two new interpolation points are chosen and the transfer function has to be evaluated only at those points. The algorithm is then extended for for time domain data based on [2].

For more details about the theory, check the paper [3].

## Main Function
```
[IdenaModel, out, wt, outequi] = adaptive_freq(H_orig, Maxite, Maxfreq, Minfreq, params)
```
***Arguments***
* `H_orig`               - TF of the original system,
* `Maxite`               - Maximum number if interations of the algorithm,
* `Maxfreq`              - Maximum frequency to be considered,
* `Minfreq`              - Minimum frequency to be considered,
* `params.tol`           - Tolerance of the Loewner framework truncation,
* `params.bw`            - Bandwidth of the filter,
* `params.NormFlag`      - A flag for norm comparison.

***Return values***

* `IdenaModel`           - Identified model using the adaptive algorithm,
* `out`                  - More information about the identified model,
* `wt`                   - Set of interpolation points used,
* `outequi`              - Model based on equidistant interpolation points.

## Examples 
The `Adaptive_freq/Examples` folder contains the two examples for the frequency domain measurments considered in the paper [3]. The `Adaptive_time/Examples` folder contains one example for the time domain measurments considered in the paper [3] based on the results in [2]. These examples reproduce the same results shown in [3].

### Penzl example 
This example is an RLC ladder circuit as part of the morwiki benchmark collection. This example has three parts:
* `Example_fom_adaptive_freq.m` This script computes a realization based on the proposed adaptive algorithm and compares it with the case where the interpolation points in an equdisitant way.
* `testNorm_fom.m` This script compares the decay in H2 error in terms of number of interpolation points for the proposed adaptive algorithm and equidistant interpolation points.
* `testbw_fom.m` This script compares the result of the algorithm for different values of the filter's bandwidth.
* `Example_fom_adaptive_freq_noisy.m` This script compares the computation of a realization using our adaptive algorithm for different levels of noise.

### Beam example 
This example is a Beam example as part of the morwiki benchmark collection. This example has two parts:
* `Example_beam_adaptive_freq.m` This script computes a realization based on the proposed adaptive algorithm and compares it with the case where the interpolation points in an equdisitant way.
* `testNorm_beam.m` This script compares the decay in H2 error in terms of number of interpolation points for the proposed adaptive algorithm and equidistant interpolation points.

### RLC example 
This example is an RLC ladder circuit [2]. This example has two parts:
* `Example_RLCSerkan_adaptive_time_limited.m` This script computes a data driven time domain realization based on the proposed adaptive algorithm.
* `Example_RLCSerkan_adaptive_time_limited_noisy.m` This script computes a data driven realization based on time domain measurments corrupted with Gaussian white noise.

## Citation
Please cite the paper [3] if you use the provided code in your research work.

## References
[1]. Mayo A. J., and Antoulas A. C., A framework for the solution of the generalized realization problem, Linear Algebra Appl. 425: 634–662, 2007.

[2]. Cherifi K., Goyal P., and Benner P., A non-intrusive method to inferring linear port-Hamiltonian realizations using time-domain data, [arXiv:2005.09371](https://arxiv.org/abs/2005.09371) 2020.

[3]. Cherifi K., Goyal P., and Benner P., Adaptive selection of interpolation points for data driven modelling, [arXiv: 2107.12950](https://arxiv.org/abs/2107.12950), 2021.
 
## Contact
Please contact [Karim Cherifi](mailto:cherifi@math.tu-berlin.de) for any queries and comments. 
___

This script has been written in MATLAB 2016b.
