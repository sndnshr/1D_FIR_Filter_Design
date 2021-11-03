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
% Examples: [h,L,err]=L2_sparse(32,1,0.5,0.6,0.004,0.00615)
%           [h,L,err]=L2_sparse(32,2,0.6,0.5,0.004,0.006)          
%           [h,L,err]=L2_sparse(64,3,[0.4 0.6], [0.265 0.73], 0.004, 0.0004);
%           [h,L,err]=L2_sparse(64,4,[0.4 0.8], [0.45 0.7], 0.01, 0.001);
function [h, L, l2_error] = L2_sparse(N1, f_type, fp_para, fa_para, mu, delta)
%% Phase I

Q = QMat(N1, f_type, fp_para, fa_para);
Ze = zeros(N1/2+1,N1/2+1);
Q1 = [Q Ze;Ze Ze];

I = eye(N1/2+1);
A = [I I;-I I];

p = pCol(N1, f_type, fp_para);
p1 = [p; mu*ones(N1/2+1,1)];

cvx_begin
    variable a(N1/2+1,1);
    variable d(N1/2+1,1);
    x =[a;d];
    minimize(x'*Q1*x + x'*p1)
    subject to
        A*x >= 0;
cvx_end

z_ind = (abs(a)<=delta);
L = sum(z_ind)*2;

%% Phase II
Q_tilde = Q;
p_tilde = p;
for i=N1/2+1:-1:1
    if z_ind(i)
    Q_tilde(:,i)=[];
    Q_tilde(i,:)=[];
    p_tilde(i,:)=[];
    end
end

a_o = -0.5*(Q_tilde\p_tilde);
af = zeros(N1/2+1,1);
k=1;
for i=1:N1/2+1
    if(~z_ind(i))
        af(i)=a_o(k);
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

figure;
plot(-1:2/511:1, fftshift(abs(fft(h,512))))
xlim([0 1])
xlabel('Normalized Frequency')
title('Amplitude response')

figure;
plot(-1:2/511:1, 20*log10(fftshift(abs(fft(h,512)))))
ylim([-50 5])
xlim([0 1])
xlabel('Normalized Frequency')
title('Amplitude response (dB)')

h_n=[a(end:-1:2)/2;a(1);a(2:end)/2];
l2_error=norm(h-h_n);
end