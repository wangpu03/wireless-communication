function [STO_est, Mag]=STO_by_correlation(y,Nfft,Ng,com_delay)
% STO estimation by maximizing the correlation between CP and rear part of OFDM symbol
% estimates STO by maximizing the correlation between CP (cyclic prefix)  
%     and rear part of OFDM symbol
% Input:  y         = Received OFDM signal including CP
%         Ng        = Number of samples in Guard Interval (CP)
%         com_delay = Common delay
% Output: STO_est   = STO estimate
%         Mag       = Correlation function trajectory varying with time

%MIMO-OFDM Wireless Communications with MATLAB¢ç   Yong Soo Cho, Jaekwon Kim, Won Young Yang and Chung G. Kang
%2010 John Wiley & Sons (Asia) Pte Ltd

N_ofdm=Nfft+Ng; 
if nargin<4, com_delay = N_ofdm/2; end
%maximum=1e-8; 
nn=0:Ng-1; yy = y(nn+com_delay)*y(nn+com_delay+Nfft)';	
maximum=abs(yy);   
for n=1:N_ofdm
   n1 = n-1;
   yy1 = y(n1+com_delay)*y(n1+com_delay+Nfft)';
   yy2 = y(n1+com_delay+Ng)*y(n1+com_delay+Nfft+Ng)';
   yy = yy-yy1+yy2;   Mag(n)=abs(yy); % Eq.(5.13)
   if (Mag(n)>maximum)
     maximum=Mag(n); STO_est = N_ofdm-com_delay-n1; 
   end
end
%if (STO_est>=Ng/2), STO_est= Ng/2-1;
% elseif (STO_est<-Ng/2), STO_est= -Ng/2;
%end    
%figure(3), plot(Mag), pause