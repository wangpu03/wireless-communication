
clc;
f1 = 900e6;				% 信号频率900MHz,信号周期为10/9 ns
N  = 10;				% 信号周期内的采样点数
Fs = N*f1;				% sampling frequency, 采样频率
T = 1/Fs;  				% sampling period, 采样周期
L = 1000*N; 			% length of signal

t = (0:L-1)*T;			% 采样时间s，fs的值越大，出来的波形失真越小
A = 4;					% 信号幅值

%% 设置高斯噪声
SNR_tag = 10;

%% 构造初始信号
source = A*sin(2*pi*f1*t);


%% 构造rayleigh信道
delay_vector = [0, 50, 110, 170, 290, 310]*1e-9; 	% Discrete delays of four-path channel (s)
gain_vector  = [0 -3.0 -10.0 -18.0 -26.0 -32.0]; 	% Average path gains (dB)			
max_Doppler_shift = 50;  					% Maximum Doppler shift of diffuse components (Hz)			
rayleigh_chan = rayleighchan(T,max_Doppler_shift,delay_vector,gain_vector);

%% 初始信号经过rayleigh信道，并保持该信道特性用于下次反射
rayleigh_chan.ResetBeforeFiltering = 1;
data_after_rayleigh = filter(rayleigh_chan,source); 

%% 计算通过rayleigh信道后的信号功率
[Pxx_hamming_after_rayleigh, F_after_rayleigh]= periodogram(data_after_rayleigh,... 
	hamming(length(data_after_rayleigh)),[],Fs,'centered', 'psd');
power_after_rayleigh = bandpower(Pxx_hamming_after_rayleigh,...
	F_after_rayleigh, 'psd');
power_after_rayleigh_db = 10*log10(power_after_rayleigh/2)



%% 初始信号经过rayleigh信道，并保持该信道特性用于下次反射
%% if reset = 0, the fading process maintains continuity from one call to the next.
rayleigh_chan.ResetBeforeFiltering = 0;
data_after_rayleigh2 = filter(rayleigh_chan,source); 

%% 计算通过rayleigh信道后的信号功率
[Pxx_hamming_after_rayleigh2, F_after_rayleigh2]= periodogram(data_after_rayleigh2,... 
	hamming(length(data_after_rayleigh2)),[],Fs,'centered', 'psd');
power_after_rayleigh2 = bandpower(Pxx_hamming_after_rayleigh2,...
	F_after_rayleigh2, 'psd');
power_after_rayleigh_db2 = 10*log10(power_after_rayleigh2/2)

isequal(data_after_rayleigh,data_after_rayleigh2)