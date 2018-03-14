%% raylei_channel.m
%%*************************************************************************************************
%% backscattering simulation
%% 测试正弦信号频谱，功率以及功率谱（三种方法）
%% 构造rayleigh衰减信道（如下表格），FFT变化，加噪声，FFT变化
%% 计算经过信道以及加噪之后信号的功率谱和信号总功率
%% @wp  2018-3-7
%%*************************************************************************************************

% Table B.2.1-2 Extended Pedestrian A model (EPA) rayleigh fading channel model
% Excess tap delay[ns/10]    	Relative power[dB]
% 0  						0.0
% 30  						-1.0
% 70						-2.0
% 90						-3.0 
% 110						-8.0
% 190						-17.2
% 410						-20.8

clc;
clear;
close all;

f1=900e6;			% 信号频率900MHz
N=20;				% 信号周期内的采样点数
Fs=N*f1;			% sampling frequency, 采样频率
T = 1/Fs;  			% sampling period, 采样周期
L = 200*N; 			% length of signal

t=(0:L-1)*T;		% 采样时间s，fs的值越大，出来的波形失真越小
A = 1;				% 信号幅值

%% 构造信号
source = A*sin(2*pi*f1*t);

figure(1);
%% display the signal
plot(t,source);
axis([0 inf -1.5 1.5]);
title('sine signal');
xlabel(['Time(f=',num2str(f1),')'])
ylabel('Amplitude/V');

data_after_fft = fft(source,L);		%对信号进行快速fourier变换
mag_fft = abs(data_after_fft);		%求得Fourier变换后的幅值

f = (1:L/2)*Fs/L;					%频率序列
%% 绘出随频率变化的振幅
figure(2);					
plot(f,mag_fft(1:L/2)*2/L);
xlabel('Frequency/Hz');
ylabel('Amplitude');title('N=10000');grid on;


%%*************************************************************************************************
%% 通过三种不同方式计算信号功率
power_theoretical = (A^2/4)*2; 	% 理论平均功率
power_theoretical_db = 10*log10(power_theoretical/2)

power_timedomain = sum(abs(source).^2)/length(source);
power_timedomain_db = 10*log10(power_timedomain/2)

% 计算该信号的功率谱
figure(3);
periodogram(source, hamming(length(source)),[],Fs,'centered', 'psd');

[Pxx_hamming, F]= periodogram(source, hamming(length(source)),[],Fs,'centered', 'psd');
power_freqdomain = bandpower(Pxx_hamming, F, 'psd');
power_freqdomain_db = 10*log10(power_freqdomain/2)
%%*************************************************************************************************

%构造rayleigh信道
delay_vector = [0, 30, 70, 90, 110, 190, 410]*1e-10; 			% Discrete delays of four-path channel (s)
gain_vector  = [0 -1.0 -2.0 -3.0 -8.0 -17.2 -20.8];  			% Average path gains (dB)
max_Doppler_shift  = 160;      									% Maximum Doppler shift of diffuse components (Hz)
chan.ResetBeforeFiltering = 0;
rayleigh_chan = rayleighchan(T,max_Doppler_shift,delay_vector,gain_vector);

%% 经过rayleigh信号
data_after_rayleigh = filter(rayleigh_chan,source); 
%%*************************************************************************************************

%% 对信道之后的进行FFT变化，并绘画出频谱图
data_after_rayleigh_fft = fft(data_after_rayleigh,L);			%对信号进行快速fourier变换
mag_after_rayleigh_fft = abs(data_after_rayleigh_fft);			%求得Fourier变换后的幅值
%% 绘出随频率变化的振幅
figure(4);					
plot(f,mag_after_rayleigh_fft(1:L/2)*2/L);
xlabel('Frequency/Hz');
ylabel('Amplitude');title('N=10000');grid on;

%% 添加噪声
SNR=10*log10(10); 
data_after_awgn=awgn(data_after_rayleigh,SNR,'measured');  
%data_after_awgn = data_after_rayleigh;


%% 对噪声之后的进行FFT变化，并绘画出频谱图
data_after_awgn_fft = fft(data_after_awgn,L);			%对信号进行快速fourier变换
mag_after_awgn_fft = abs(data_after_awgn_fft);			%求得Fourier变换后的幅值
%% 绘出随频率变化的振幅
figure(5);					
plot(f,mag_after_awgn_fft(1:L/2)*2/L);
xlabel('Frequency/Hz');
ylabel('Amplitude');title('N=10000');grid on;

%%*************************************************************************************************

% 在计算的过程中已经画出图形，后面都是对图形的调整
[Pxx_hamming_after_awgn, F_after_awgn]= periodogram(data_after_awgn, hamming(length(data_after_awgn)),[],Fs,'centered', 'psd');
power_freqdomain_after_awgn = bandpower(Pxx_hamming_after_awgn, F_after_awgn, 'psd');
power_freqdomain_after_awgn_db = 10*log10(power_freqdomain_after_awgn/2)

figure;
plot(t,real(data_after_awgn));
%axis([0 inf -1.5 1.5]);
title('sine signal after channel');
xlabel(['Time(f=',num2str(f1),')'])
ylabel('Amplitude/V');