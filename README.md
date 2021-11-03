# 1D FIR Sparse Filter Design Using Optimization

This work implements the 1D FIR sparse filter design technique proposed in the following paper for lowpass, highpass, bandpass and bandstop filters.

\[1\] W.-S. Lu and T. Hinamoto, "Digital flters with sparse coefficients," in Proceedings of 2010 IEEE International Symposium on Circuits and Systems. IEEE, 2010, pp. 169-172.

## Functions

- L2_sparse(N1, f_type, fp_para, fa_para, mu, delta) 

Designs the desired standard filters in the least squares sense. Outputs the impulse response, number of zeros, and the L-2 error.

- L_inf_sparse(N1, f_type, fp_para, fa_para, mu, delta) 

Designes the desired standard filters in the minimax sense. Outputs the impulse response, number of zeros, and the L_inf error

- QMat(N1, f_type, fp_para, fa_para)

Supporting function for L2_sparse. Calculates the Hessian matrix.

- pCol(N1, f_type, fp_para)

Supporting function for L2_sparse. Calculates the gradient vector needed for the first phase of L2 sparse filter design.

- freqGrids(N1, f_type, fp_para, fa_para)

Supporting function for L_inf_sparse. Outputs a frequency grid of passband and stopband frequencies.

### Function parameters

- N1: filter order (even integer).
- f_type: 1 for lowpass, 2 for highpass, 3 for bandpass, and 4 for bandstop.
- fp_para: Normalized passband frequency fp for types 1 and 2, [fp1 fp2] for types 3 and 4.
- fa_para: Normalized stopband frequency fa for types 1 and 2, [fa1 fa2] for types 3 and 4.
- mu: regularization parameter.
- delta: threshold for hardthresholding.