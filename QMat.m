% Supporting function for L2_sparse.m
% Inputs: N1 -- filter order
%         f_type -- 1 for lowpass, 2 for highpass, 3 for bandpass,
%                   and 4 for bandstop.
%         fp_para -- Normalized passband frequency fp for types 1 and 2,
%                   [fp1 fp2] for types 3 and 4
%         fa_para -- Normalized stopband frequency fa for types 1 and 2,
%                   [fa1 fa2] for types 3 and 4
% Output: Q -- Hessian matrix for the first phase of least-squares sparse filter design
function Q = QMat(N1, f_type, fp_para, fa_para)
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

Q = zeros(N1/2+1,N1/2+1);
if f_type == 1
    for i=1:N1/2+1
        for j=1:N1/2+1
            if i==1 && j==1
                Q(1,1)=2*(pi+wp-wa);
            elseif i==j
                Q(i,i)= pi + wp - wa + sin(2*(i-1)*wp)/(2*(i-1)) - sin(2*(i-1)*wa)/(2*(i-1));
            else
                Q(i,j)= sin((i+j-2)*wp)/(i+j-2) - sin((i+j-2)*wa)/(i+j-2) + sin((i-j)*wp)/(i-j) - sin((i-j)*wa)/(i-j);
            end
        end
    end
elseif f_type == 2
    for i=1:N1/2+1
        for j=1:N1/2+1
            if i==1 && j==1
                Q(1,1)=2*(pi+wa-wp);
            elseif i==j
                Q(i,i)= pi + wa - wp + sin(2*(i-1)*wa)/(2*(i-1)) - sin(2*(i-1)*wp)/(2*(i-1));
            else
                Q(i,j)= sin((i+j-2)*wa)/(i+j-2) - sin((i+j-2)*wp)/(i+j-2) + sin((i-j)*wa)/(i-j) - sin((i-j)*wp)/(i-j);
            end
        end
    end
elseif f_type == 3
    for i=1:N1/2+1
        for j=1:N1/2+1
            if i==1 && j==1
                Q(1,1)=2*(pi + wp2 - wp1 + wa1 - wa2);
            elseif i==j
                Q(i,i)= pi + wp2 - wp1 + wa1 - wa2 + sin(2*(i-1)*wp2)/(2*(i-1)) - sin(2*(i-1)*wp1)/(2*(i-1))...
                    + sin(2*(i-1)*wa1)/(2*(i-1)) - sin(2*(i-1)*wa2)/(2*(i-1));
            else
                Q(i,j)= sin((i+j-2)*wp2)/(i+j-2) - sin((i+j-2)*wp1)/(i+j-2) + sin((i-j)*wp2)/(i-j) - sin((i-j)*wp1)/(i-j)...
                    + sin((i+j-2)*wa1)/(i+j-2) - sin((i+j-2)*wa2)/(i+j-2) + sin((i-j)*wa1)/(i-j) - sin((i-j)*wa2)/(i-j);
            end
        end
    end
else
    for i=1:N1/2+1
        for j=1:N1/2+1
            if i==1 && j==1
                Q(1,1)=2*(pi + wp1 - wp2 + wa2 - wa1);
            elseif i==j
                Q(i,i)= pi + wp1 - wp2 + wa2 - wa1 + sin(2*(i-1)*wp1)/(2*(i-1)) - sin(2*(i-1)*wp2)/(2*(i-1))...
                    + sin(2*(i-1)*wa2)/(2*(i-1)) - sin(2*(i-1)*wa1)/(2*(i-1));
            else
                Q(i,j)= sin((i+j-2)*wp1)/(i+j-2) - sin((i+j-2)*wp2)/(i+j-2) + sin((i-j)*wp1)/(i-j) - sin((i-j)*wp2)/(i-j)...
                    + sin((i+j-2)*wa2)/(i+j-2) - sin((i+j-2)*wa1)/(i+j-2) + sin((i-j)*wa2)/(i-j) - sin((i-j)*wa1)/(i-j);
            end
        end
    end
end
end