function CFO_est=CFO_CP(y,Nfft,Ng)
% Time-domain CFO estimation based on CP (Cyclic Prefix)

%MIMO-OFDM Wireless Communications with MATLAB¢ç   Yong Soo Cho, Jaekwon Kim, Won Young Yang and Chung G. Kang
%2010 John Wiley & Sons (Asia) Pte Ltd

nn=1:Ng; CFO_est = angle(y(nn+Nfft)*y(nn)')/(2*pi);  % Eq.(5.27)  
% MSE = MSE + (CFO_est-CFO)*(CFO_est-CFO);