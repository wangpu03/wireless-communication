%  例2：x=0.5*sin(2*pi*15*t)+2*sin(2*pi*40*t),fs=100Hz,绘制：
% （1）数据个数N=32，FFT所用的采样点数NFFT=32；
% （2）N=32，NFFT=128；
% （3）N=136，NFFT=128；
% （4）N=136，NFFT=512。


%  结论：
% （1）当数据个数和FFT采用的数据个数均为32时，频率分辨率较低，但没有由于
%  添零而导致的其他频率成分。
% （2）由于在时间域内信号加零，致使振幅谱中出现很多其他成分，这是加零造成的。
%  其振幅由于加了多个零而明显减小。
% （3）FFT程序将数据截断，这时分辨率较高。
% （4）也是在数据的末尾补零，但由于含有信号的数据个数足够多，FFT振幅谱也基本不受影响。
% 
%  对信号进行频谱分析时，数据样本应有足够的长度，一般FFT程序中所用数据点数与原含有
%  信号数据点数相同，这样的频谱图具有较高的质量，可减小因补零或截断而产生的影响。


clc;
clear;
fs = 100;
Ndata = 32;
N = 32; 		%FFT的长度
n = 0:Ndata-1;
t = n/fs;  		%数据对应的时间序列

x = 0.5*sin(2*pi*15*t)+2*sin(2*pi*40*t);

y = fft(x,N);
mag = abs(y);
f = (0:N-1)*fs/N;

subplot(2,2,1);
plot(f(1:N/2),mag(1:N/2)*2/N);	%绘出Nyquist频率之前的振幅
xlabel('F/Hz');ylabel('Amplitude');
title('Ndata=32 Nfft=32');grid on;

Ndata = 32;
N = 128; 		%FFT的长度
n = 0:Ndata-1;
t = n/fs;  		%数据对应的时间序列

x = 0.5*sin(2*pi*15*t)+2*sin(2*pi*40*t);

y = fft(x,N);
mag = abs(y);
f = (0:N-1)*fs/N;

subplot(2,2,2);
plot(f(1:N/2),mag(1:N/2)*2/N);	%绘出Nyquist频率之前的振幅
xlabel('F/Hz');ylabel('Amplitude');
title('Ndata=32 Nfft=128');grid on;

Ndata = 136;
N = 128; 		%FFT的长度
n = 0:Ndata-1;
t = n/fs;  		%数据对应的时间序列

x = 0.5*sin(2*pi*15*t)+2*sin(2*pi*40*t);

y = fft(x,N);
mag = abs(y);
f = (0:N-1)*fs/N;

subplot(2,2,3);
plot(f(1:N/2),mag(1:N/2)*2/N);	%绘出Nyquist频率之前的振幅
xlabel('F/Hz');ylabel('Amplitude');
title('Ndata=136 Nfft=128');grid on;


Ndata = 136;
N = 512; 		%FFT的长度
n = 0:Ndata-1;
t = n/fs;  		%数据对应的时间序列

x = 0.5*sin(2*pi*15*t)+2*sin(2*pi*40*t);

y = fft(x,N);
mag = abs(y);
f = (0:N-1)*fs/N;

subplot(2,2,4);
plot(f(1:N/2),mag(1:N/2)*2/N);	%绘出Nyquist频率之前的振幅
xlabel('F/Hz');ylabel('Amplitude');
title('Ndata=136 Nfft=512');grid on;
