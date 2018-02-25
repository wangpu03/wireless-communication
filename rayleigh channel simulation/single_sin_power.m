% measuring the power of deterministic periods signals

clc;
close all;
Fs = 1024;   					% 采样频率
t  = 0:1/Fs:1-(1/Fs);			% 时间序列
A  = 1;   						% Vpeak
F1 = 128; 						% Hz，信号频率
x  = A*sin(2*pi*t*F1);

idx = 1:128;
% figure;
% plot(t(idx),x(idx));
% grid;
% ylabel('Amplitude');
% xlabel('Time (sec)');
% axis tight;

power_theoretical = (A^2/4)*2 	% 理论平均功率
power_theoretical_db = 10*log10(power_theoretical/2)

%% Measuring the Power of a Single Sinusoid
% figure;
% % 在计算的过程中已经画出图形，后面都是对图形的调整
% periodogram(x, hamming(length(x)),[],Fs,'centered', 'power');
% v = axis;
% axis([v(1) v(2) -10 -5.5])
% hgcf = gcf;						% 图形句柄
% hgcf.Color = [1 1 1];


% %% Estimating the Power of a Single Sinusoid via PSD
figure
power_freqdomain = periodogram(x, hamming(length(x)), [], Fs, 'centered', 'psd');
v = axis;
axis([v(1) v(2) -10 -5.5])
hgcf = gcf;
hgcf.Color = [1 1 1];

%% we notice in thsi plot is that teh peaks of the spectrum plot do not
%% have the same height as when we plotted the power spectrum.
%% because when taking Power Spectral density (PSD) measurement it's the erea 
%% under the curve that matters. 
%% 即峰值线面的曲线部分也会占据一部分能量，所以导致峰值下降。


[Pxx_hamming, F] = periodogram(x, hamming(length(x)), [], Fs, 'psd');
%% 经过功率谱变换后，可以利用这个公式计算信号频率域上的功率
power_freqdomain = bandpower(Pxx_hamming, F, 'psd');

%% 下面公式是计算信号时间域上的公式
power_timedomain = sum(abs(x).^2)/length(x);



