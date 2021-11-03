% Supporting function for L2_sparse.m
% Inputs: N1 -- filter order
%         f_type -- 1 for lowpass, 2 for highpass, 3 for bandpass,
%                   and 4 for bandstop.
%         fp_para -- Normalized passband frequency fp for types 1 and 2,
%                   [fp1 fp2] for types 3 and 4
% Output: p -- Gradient vector for the first phase of least-squares sparse
%              filter design
function p = pCol(N1, f_type, fp_para)

if f_type == 1 % lowpass
   wp = pi*fp_para;
   p = -4*([wp sin((1:1:N1/2)*wp)./((1:1:N1/2))]');
elseif f_type == 2 % highpass
   wp = pi*fp_para;
   p = 4*([(wp-pi) sin((1:1:N1/2)*wp)./((1:1:N1/2))]');
elseif f_type == 3 % bandpass
   wp1 = pi*fp_para(1);
   wp2 = pi*fp_para(2);
   p = -4*([(wp2-wp1) (sin((1:1:N1/2)*wp2)-sin((1:1:N1/2)*wp1))./((1:1:N1/2))]');
elseif f_type == 4 % bandstop
   wp1 = pi*fp_para(1);
   wp2 = pi*fp_para(2);
   p = -4*([(pi + wp1-wp2) (sin((1:1:N1/2)*wp1)-sin((1:1:N1/2)*wp2))./((1:1:N1/2))]');
end
end