% Inputs: N1 -- filter order
%         f_type -- 1 for lowpass, 2 for highpass, 3 for bandpass,
%                   and 4 for bandstop.
%         fp_para -- Normalized passband frequency fp for types 1 and 2,
%                   [fp1 fp2] for types 3 and 4
%         fa_para -- Normalized stopband frequency fa for types 1 and 2,
%                   [fa1 fa2] for types 3 and 4
%         mu -- Regularization parameter
%         delta -- Threshold for hard thresholding
% Output: h -- filter's impulse response (coefficients).
% Examples: [h,L,err] = L_inf_sparse(64,1,0.5,0.6,0.07,0.003);
%           [h,L,err] = L_inf_sparse(64,2,0.6,0.5,0.07,0.003);        
%           [h,L,err] = L_inf_sparse(64,3,[0.4 0.6],[0.265 0.73],0.004,0.0004);
%           [h,L,err] = L_inf_sparse(64,4,[0.4 0.8],[0.5 0.7],0.004,0.003);
function [h, L, l_inf_error] = L_inf_sparse(N1, f_type, fp_para, fa_para, mu, delta)

[Ap, As] = freq_grids(N1, f_type, fp_para, fa_para);

cvx_begin
    variable n;
    variable u(N1/2+1,1);
    variable v(N1/2+1,1);
    x = [n; u; v];
    c = [1; ones(N1/2+1,1)*mu; ones(N1/2+1,1)*mu];
    
    minimize(c'*x)
    subject to:
        abs(Ap*(u - v) - 1) <= n;
        abs(As*(u - v)) <= n;
        u >= 0;
        v >= 0;
cvx_end

a = u-v;

%% Phase II

z_ind = (abs(a)<=delta);
L = sum(z_ind)*2;
for i=N1/2+1:-1:1
    if z_ind(i)
    Ap(:,i)=[];
    As(:,i)=[];
    end
end

A = [Ap; As];
b = [ones(size(Ap,1),1); zeros(size(As,1),1)];
cvx_begin
    variable n;
    variable a1(N1/2+1-L/2,1);
    minimize(n)
    subject to:
    abs(A*a1 - b) <= n;
cvx_end

af = zeros(N1/2+1,1);
k=1;
for i=1:N1/2+1
    if(~z_ind(i))
        af(i)=a1(k);
        k=k+1;
    end
end

%impulse response h
h = zeros(N1+1,1);
h(N1/2+1)=af(1);
for i=N1/2:-1:1
   h(i)=0.5*af(N1/2+2-i); 
   h(N1+2-i)=h(i);
end

figure;
stem(-N1/2:1:N1/2,h)
title('Impulse response');

figure
plot(-1:2/511:1, fftshift(abs(fft(h,512))))
xlim([0 1])
xlabel('Normalized Frequency')
title('Amplitude response')

figure;
plot(-1:2/511:1, 20*log10(fftshift(abs(fft(h,512)))))
ylim([-60 5])
xlim([0 1])
xlabel('Normalized Frequency')
title('Amplitude response (dB)')

h_n=[a(end:-1:2)/2;a(1);a(2:end)/2];
l_inf_error=max(h-h_n);
end