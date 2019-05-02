function xp=add_pilot(x,Nfft,Nps)
% Nps : Pilot spacing

%MIMO-OFDM Wireless Communications with MATLAB¢ç   Yong Soo Cho, Jaekwon Kim, Won Young Yang and Chung G. Kang
%2010 John Wiley & Sons (Asia) Pte Ltd

if nargin<3,  Nps=4;  end
Np=Nfft/Nps;  xp=x; % Number of pilots and an OFDM signal including pilot signal
for k=1:Np
   xp((k-1)*Nps+1)= exp(j*pi*(k-1)^2/Np);  % Pilot boosting with an even Np
   %xp((k-1)*Nps+1)= exp(j*pi*(k-1)*k/Np);  % Pilot boosting with an odd Np
end
% CAZAC (Constant Amplitude Zero AutoCorrelation) sequence --> pilot