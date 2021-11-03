% Supporting function for L_inf_sparse.m
% Inputs: N1 -- filter order
%         f_type -- 1 for lowpass, 2 for highpass, 3 for bandpass,
%                   and 4 for bandstop.
%         fp_para -- Normalized passband frequency fp for types 1 and 2,
%                   [fp1 fp2] for types 3 and 4
%         fa_para -- Normalized stopband frequency fa for types 1 and 2,
%                   [fa1 fa2] for types 3 and 4
% Output: Ap -- Passband frequency grid
%         As -- Stopband frequency grid
function [Ap, As] = freq_grids(N1, f_type, fp_para, fa_para)

if f_type == 1 % lowpass
   wp = pi*fp_para;
   wa = pi*fa_para;
elseif f_type == 2 % highpass
   wp = pi*fp_para;
   wa = pi*fa_para;
elseif f_type == 3 % bandpass
   wp1 = pi*fp_para(1);
   wp2 = pi*fp_para(2);
   wa1 = pi*fa_para(1);
   wa2 = pi*fa_para(2);
elseif f_type == 4 % bandstop
   wp1 = pi*fp_para(1);
   wp2 = pi*fp_para(2);
   wa1 = pi*fa_para(1);
   wa2 = pi*fa_para(2);
end

P = N1*15;
nn = 1:N1/2;

if f_type == 1
    ww1 = 0:wp/(P-1):wp;
    [n1,w1] = meshgrid(nn,ww1);
    Ap = [ones(P,1) , cos(w1.*n1)];

    ww2 = wa:(pi-wa)/(P-1):pi;
    [n2,w2] = meshgrid(nn,ww2);
    As = [ones(P,1) , cos(w2.*n2)];

elseif f_type == 2   
    ww1 = wp:(pi-wp)/(P-1):pi;
    [n1,w1] = meshgrid(nn,ww1);
    Ap = [ones(P,1) , cos(w1.*n1)];

    ww2 = 0:wa/(P-1):wa;
    [n2,w2] = meshgrid(nn,ww2);
    As = [ones(P,1) , cos(w2.*n2)];
elseif f_type == 3
    ww1 = 0:wa1/(P-1):wa1;
    [n1,w1] = meshgrid(nn,ww1);
    As1 = [ones(P,1) , cos(w1.*n1)];

    ww2 = wp1:(wp2-wp1)/(P-1):wp2;
    [n2,w2] = meshgrid(nn,ww2);
    Ap = [ones(P,1) , cos(w2.*n2)];

    ww3 = wa2:(pi-wa2)/(P-1):pi;
    [n3,w3] = meshgrid(nn,ww3);
    As2 = [ones(P,1) , cos(w3.*n3)];
    
    As = [As1; As2];
else
    ww1 = 0:wp1/(P-1):wp1;
    [n1,w1] = meshgrid(nn,ww1);
    Ap1 = [ones(P,1) , cos(w1.*n1)];

    ww2 = wa1:(wa2-wa1)/(P-1):wa2;
    [n2,w2] = meshgrid(nn,ww2);
    As = [ones(P,1) , cos(w2.*n2)];

    ww3 = wp2:(pi-wp2)/(P-1):pi;
    [n3,w3] = meshgrid(nn,ww3);
    Ap2 = [ones(P,1) , cos(w3.*n3)];
    
    Ap = [Ap1; Ap2];
end
end