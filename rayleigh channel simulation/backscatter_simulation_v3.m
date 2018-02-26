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
L = 10000*N; 			% length of signal

t = (0:L-1)*T;			% 采样时间s，fs的值越大，出来的波形失真越小
A = 4;					% 信号幅值

%% 设置高斯噪声
SNR_tag = 10;

%% 构造初始信号
source = A*sin(2*pi*f1*t);

%% 计算初始信号的功率
[Pxx_hamming, F]= periodogram(source, hamming(length(source)),[],Fs,'centered', 'psd');
power_source = bandpower(Pxx_hamming, F, 'psd');
power_source_db = 10*log10(power_source/2)


%% 构造rayleigh信道
delay_vector = [0, 50, 110, 170, 290, 310]*1e-9; 	% Discrete delays of four-path channel (s)
gain_vector  = [0 -3.0 -10.0 -18.0 -26.0 -32.0]; 	% Average path gains (dB)			
max_Doppler_shift = 50;  					% Maximum Doppler shift of diffuse components (Hz)			
rayleigh_chan = rayleighchan(T,max_Doppler_shift,delay_vector,gain_vector);

%% 初始信号经过rayleigh信道，并保持该信道特性用于下次反射
rayleigh_chan.ResetBeforeFiltering = 0;
data_after_rayleigh = filter(rayleigh_chan,source); 

%% 计算通过rayleigh信道后的信号功率
[Pxx_hamming_after_rayleigh, F_after_rayleigh]= periodogram(data_after_rayleigh,... 
	hamming(length(data_after_rayleigh)),[],Fs,'centered', 'psd');
power_after_rayleigh = bandpower(Pxx_hamming_after_rayleigh,...
	F_after_rayleigh, 'psd');
power_after_rayleigh_db = 10*log10(power_after_rayleigh/2)


%% 添加高斯噪声
data_tag = awgn(data_after_rayleigh,SNR_tag,'measured');

%% 计算在tag端接收信号的功率
[Pxx_hamming_tag, F_tag] = periodogram(data_tag,hamming(length(data_tag)),[],Fs,'centered','psd');
power_tag = bandpower(Pxx_hamming_tag,F_tag,'psd');
power_tag_db = 10*log10(power_tag/2)

%%************************************************************
%% 反射路径 
%% 发射因子为0.5
coeffi = 1.0;
data_back = data_tag.*coeffi;

%% 反射后，信号经过再次经过rayleigh信道
rayleigh_chan.ResetBeforeFiltering = 0;
back_after_rayleigh = filter(rayleigh_chan,data_back); 

%% 计算反射之后经过rayleigh信道后的信号功率
[Pxx_hamming_after_back, F_after_back]= periodogram(back_after_rayleigh,...
	hamming(length(back_after_rayleigh)),[],Fs,'centered', 'psd');
power_after_back = bandpower(Pxx_hamming_after_back, F_after_back, 'psd');
power_after_back_db = 10*log10(power_after_back/2)


%% 添加高斯噪声
data_pt = awgn(back_after_rayleigh,SNR_tag,'measured');

%% 计算在tag端接收信号的功率
[Pxx_hamming_pt, F_pt] = periodogram(data_pt,hamming(length(data_pt)),[],Fs,'centered','psd');
power_pt = bandpower(Pxx_hamming_pt,F_pt,'psd');
power_pt_db = 10*log10(power_pt/2)



%% 绘制最终PT端接收的信号图形
figure;
plot(t(1:1000),real(data_pt(1:1000)));
% axis([0 inf -1.5 1.5]);
title('sine signal after channel');
xlabel(['Time(f=',num2str(f1),')'])
ylabel('Amplitude/V');

