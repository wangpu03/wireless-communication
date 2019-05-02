function y_CFO=add_CFO(y,CFO,Nfft)
% To add an arbitrary frequency offset
% Input: y    = Time-domain received signal
%        dCFO = FFO (fractional CFO) + IFO (integral CFO)
%        Nfft = FFT size;

%MIMO-OFDM Wireless Communications with MATLAB¢ç   Yong Soo Cho, Jaekwon Kim, Won Young Yang and Chung G. Kang
%2010 John Wiley & Sons (Asia) Pte Ltd

nn=0:length(y)-1; y_CFO = y.*exp(j*2*pi*CFO*nn/Nfft);
% plot(real(y_CFO),imag(y_CFO),'.')
    