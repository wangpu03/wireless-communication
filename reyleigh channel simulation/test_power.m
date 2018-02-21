clear all
Fs=1000;
n=0:1/Fs:1;
x=cos(2*pi*40*n)+3*cos(2*pi*100*n)+randn(size(n));

nfft=1024;
window=boxcar(length(n));
[Pxx,f]=periodogram(x,window,nfft,Fs);
P=10*log10(Pxx);
plot(f,P);
hold on;
Pxx_1=abs(fft(x,nfft)).^2/length(n);
t=0:round(nfft/2-1);
f=t*Fs/nfft;
P_1=10*log10(Pxx_1(t+1));
plot(f,P_1,'r');
legend('periodogram','??');
title('?????????');