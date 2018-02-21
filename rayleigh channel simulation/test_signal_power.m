clear;
Fs=1000; %采样频率
n=0:1/Fs:1;
%产生含有噪声的序列
xn=cos(2*pi*40*n)+3*cos(2*pi*100*n)+randn(size(n));
window=boxcar(length(xn)); %矩形窗
nfft=1024;
[Pxx,f]=periodogram(xn,window,nfft,Fs); %直接法
subplot(2,2,1);
plot(f,10*log10(Pxx));

clear;
Fs=1000; %采样频率
n=0:1/Fs:1;
%产生含有噪声的序列
xn=cos(2*pi*40*n)+3*cos(2*pi*100*n)+randn(size(n));
nfft=1024;
cxn=xcorr(xn,'unbiased'); %计算序列的自相关函数
CXk=fft(cxn,nfft);
Pxx=abs(CXk);
index=0:round(nfft/2-1);
k=index*Fs/nfft;
plot_Pxx=10*log10(Pxx(index+1));
subplot(2,2,2);
plot(k,plot_Pxx);

clear;
Fs=1000;
n=0:1/Fs:1;
xn=cos(2*pi*40*n)+3*cos(2*pi*100*n)+randn(size(n));
nfft=1024;
window=boxcar(100); %矩形窗
window1=hamming(100); %海明窗
window2=blackman(100); %blackman窗
noverlap=20; %数据无重叠
range='half'; %频率间隔为[0 Fs/2]，只计算一半的频率
[Pxx,f]=pwelch(xn,window,noverlap,nfft,Fs,range);
[Pxx1,f]=pwelch(xn,window1,noverlap,nfft,Fs,range);
[Pxx2,f]=pwelch(xn,window2,noverlap,nfft,Fs,range);
plot_Pxx=10*log10(Pxx);
plot_Pxx1=10*log10(Pxx1);
plot_Pxx2=10*log10(Pxx2);
subplot(2,2,3);
plot(f,plot_Pxx);
subplot(2,2,4);
plot(f,plot_Pxx1);
figure(2)
plot(f,plot_Pxx2); 
clc;