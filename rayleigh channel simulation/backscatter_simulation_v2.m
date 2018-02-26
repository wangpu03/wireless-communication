% backscattering simulation
% Table B.2.1-2 Extended Pedestrian A model (EPA) rayleigh fading channel model
% Excess tap delay[ns/10]    	Relative power[dB]
% 0  						0.0
% 50  						-3.0
% 110						-10.0
% 170						-18.0 
% 290						-26.0
% 310						-32.0


%% 该部分测试正弦信号经过rayleigh衰减信道，和经过rayleigh衰减信道后并添加部分噪声后的  
%% 信号功率的比较

clc;
clear;
close all;

f1 = 900e6;				% 信号频率900MHz,信号周期为10/9 ns
N  = 10;				% 信号周期内的采样点数
Fs = N*f1;				% sampling frequency, 采样频率
T = 1/Fs;  				% sampling period, 采样周期
L = 20000*N; 				% length of signal

t = (0:L-1)*T;			% 采样时间s，fs的值越大，出来的波形失真越小
A = 1;					% 信号幅值

%% 构造信号
source = A*sin(2*pi*f1*t);
[Pxx_hamming, F]= periodogram(source, hamming(length(source)),[],Fs,'centered', 'psd');
power_freqdomain = bandpower(Pxx_hamming, F, 'psd');
power_freqdomain_db = 10*log10(power_freqdomain/2)


%% 构造rayleigh信道
% Discrete delays of four-path channel (s)
delay_vector = [0, 50, 110, 170, 290, 310]*1e-10; 

%% Average path gains (dB)			
gain_vector  = [0 -3.0 -10.0 -18.0 -26.0 -32.0]; 

%% Maximum Doppler shift of diffuse components (Hz)			
max_Doppler_shift = 50;      

rayleigh_chan = rayleighchan(T,max_Doppler_shift,delay_vector,gain_vector);

%% 经过rayleigh信道，且将保持该信道特性
rayleigh_chan.ResetBeforeFiltering = 0;
data_after_rayleigh = filter(rayleigh_chan,source); 


%% 功率谱，并计算整个信号的功率
[Pxx_hamming_after_rayleigh, F_after_rayleigh]= periodogram(data_after_rayleigh,... 
	hamming(length(data_after_rayleigh)),[],Fs,'centered', 'psd');
power_freqdomain_after_rayleigh = bandpower(Pxx_hamming_after_rayleigh,...
	F_after_rayleigh, 'psd');
power_freqdomain_after_rayleigh_db = 10*log10(power_freqdomain_after_rayleigh/2)


%% 添加高斯噪声，
SNR = 10*log10(100); 
data_after_awgn=awgn(data_after_rayleigh,SNR,'measured');

%% 功率谱，并计算tag出接收到的整个信号功率
[Pxx_hamming_after_awgn, F_after_awgn]= periodogram(data_after_awgn,...
	hamming(length(data_after_awgn)),[],Fs,'centered', 'psd');
power_freqdomain_after_awgn = bandpower(Pxx_hamming_after_awgn, F_after_awgn, 'psd');
power_freqdomain_after_awgn_db = 10*log10(power_freqdomain_after_awgn/2)


%%************************************************************
%% 反射路径 

%% 发射因子为0.5
coeffi = 0.5;
data_backscatter = data_after_awgn.*coeffi;

rayleigh_chan.ResetBeforeFiltering = 0;
backscatter_after_rayleigh = filter(rayleigh_chan,data_backscatter); 

%% 功率谱，并计算backscatter回去的整个信号功率
[Pxx_hamming_after_back, F_after_back]= periodogram(data_backscatter,...
	hamming(length(data_backscatter)),[],Fs,'centered', 'psd');
power_freqdomain_after_back = bandpower(Pxx_hamming_after_back, F_after_back, 'psd');
power_freqdomain_after_back_db = 10*log10(power_freqdomain_after_back/2)

%% 返回去的信道中添加高斯噪声，
SNR_back = 10*log10(100); 
back_after_awgn=awgn(backscatter_after_rayleigh,SNR_back,'measured');

%% 功率谱，并计算接收方的信号功率
[Pxx_hamming_after_back_agwn, F_after_back_agwn]= periodogram(back_after_awgn,...
	hamming(length(back_after_awgn)),[],Fs,'centered', 'psd');
power_freqdomain_after_back_agwn = bandpower(Pxx_hamming_after_back_agwn,...
	F_after_back_agwn, 'psd');
power_freqdomain_after_back_agwn_db = 10*log10(power_freqdomain_after_back_agwn/2)


figure;
plot(t(1:2000),real(back_after_awgn(1:2000)));
%axis([0 inf -1.5 1.5]);
title('sine signal after channel');
xlabel(['Time(f=',num2str(f1),')'])
ylabel('Amplitude/V');

