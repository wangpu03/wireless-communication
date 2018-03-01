% Table is a rayleigh fading channel model in indoor
% Excess tap delay[ns]    	Relative power[dB]
% 0  						0.0
% 50  						-3.0
% 110						-10.0
% 170						-18.0 
% 290						-26.0
% 310						-32.0
%% 之后的信号都这个模型的基础上进行的修改

%% 仿真方案二，设备之间的密钥建立，但是在tag能够区分整个信号时，即不叠加两个信号

clc;
clear;
close all;

f1 = 900e6;				% 信号频率900MHz,信号周期为10/9 ns
N  = 10;				% 信号周期内的采样点数
Fs = N*f1;				% sampling frequency, 采样频率
T = 1/Fs;  				% sampling period, 采样周期
L = 10000*N; 			% length of signal

t = (0:(L-1))*T;			% 采样时间s，fs的值越大，出来的波形失真越小
A = 4;					% 信号幅值



%% 构造调制1和-1序列
f_seq = 22.5e6;   %%	45MHz
value = 1;
for index = 1:L
		base_signal(index) = value;
	if mod(index, (Fs/f_seq)/2) == 0
		value = 0-value;
	end
end

% base_signal=square(2*pi*f_seq*t);
% sum(base_signal)



%% 构造初始信号
source = A*sin(2*pi*f1*t);

plot(t(1:200),source(1:200),t(1:200),base_signal(1:200))

%% 计算初始信号的功率
[Pxx_hamming, F]= periodogram(source, hamming(length(source)),[],Fs,'centered', 'psd');
power_source = bandpower(Pxx_hamming, F, 'psd');
power_source_db = 10*log10(power_source/2);


num = 1:1000;

power_source_array(num) = power_source_db;

power_tag_array(num) = 0;

power_pt_array(num) = 0;

for index = num

	%% 测试PT传输信号到Alice时的信号功率
	%% 构造rayleigh信道一，表示PT与Alice之间的信道
	delay_vector_a = [0, 50, 110, 170, 290, 310]*1e-9; 	% Discrete delays of four-path channel (s)
	gain_vector_a  = [0 -3.0 -10.0 -18.0 -26.0 -32.0]; 	% Average path gains (dB)			
	max_Doppler_shift_a = 50;  					% Maximum Doppler shift of diffuse components (Hz)			
	rayleigh_chan_a = rayleighchan(T,max_Doppler_shift_a,delay_vector_a,gain_vector_a);

	%% 初始信号经过rayleigh信道，并保持该信道特性用于下次反射
	rayleigh_chan_a.ResetBeforeFiltering = 0;
	data_after_rayleigh_a = filter(rayleigh_chan_a,source); 

	%% 设置高斯噪声与信号的信噪比
	SNR_tag_a = 8;
	%% 添加高斯噪声
	data_tag_a = awgn(data_after_rayleigh_a,SNR_tag_a,'measured');

	%% 计算在tag端接收信号的功率
	[Pxx_hamming_tag_a, F_tag_a] = periodogram(data_tag_a,hamming(length(data_tag_a)),[],Fs,'centered','psd');
	power_tag_a = bandpower(Pxx_hamming_tag_a,F_tag_a,'psd');
	power_tag_a_db = 10*log10(power_tag_a/2);

	power_tag_a_array(index) = power_tag_a;


	%%**********************************************************************************************
	%% 测试PT传输信号到Bob时的信号功率
	%% 构造rayleigh信道二，表示PT与Bob之间的信道
	delay_vector_b = [0, 55, 120, 190, 310, 330]*1e-9; 	% Discrete delays of four-path channel (s)
	gain_vector_b  = [0 -4.0 -11.0 -20.0 -28.0 -36.0]; 	% Average path gains (dB)			
	max_Doppler_shift_b = 55;  					% Maximum Doppler shift of diffuse components (Hz)			
	rayleigh_chan_b = rayleighchan(T,max_Doppler_shift_b,delay_vector_b,gain_vector_b);

	%% 初始信号经过rayleigh信道，并保持该信道特性用于下次反射
	rayleigh_chan_b.ResetBeforeFiltering = 0;
	data_after_rayleigh_b = filter(rayleigh_chan_b,source); 

	%% 设置高斯噪声与信号的信噪比
	SNR_tag_b = 9;
	%% 添加高斯噪声
	data_tag_b = awgn(data_after_rayleigh_b,SNR_tag_b,'measured');

	%% 计算在tag端接收信号的功率
	[Pxx_hamming_tag_b, F_tag_b] = periodogram(data_tag_b,hamming(length(data_tag_b)),[],Fs,'centered','psd');
	power_tag_b = bandpower(Pxx_hamming_tag_b,F_tag_b,'psd');
	power_tag_b_db = 10*log10(power_tag_b/2);

	power_tag_b_array(index) = power_tag_b;


	%%**********************************************************************************************
	%% Alice的反射路径，将信号反射给Bob
	coeffi_a = 0.9;						%% 反射因子
	data_a_coeffi = data_tag_a.*coeffi_a;

	%% 进行ASK调制，使用已经设置好的方波序列
	data_a_back = data_a_coeffi.*base_signal;

	%% 构造rayleigh信道，表示Alice与Bob之间的信道
	delay_vector_ab = [0, 4, 10, 19, 31, 43]*1e-9; 	% Discrete delays of four-path channel (s)
	gain_vector_ab  = [0 -1.0 -5.0 -10.0 -18.0 -26.0]; 	% Average path gains (dB)			
	max_Doppler_shift_ab = 100;  					% Maximum Doppler shift of diffuse components (Hz)			
	rayleigh_chan_ab = rayleighchan(T,max_Doppler_shift_ab,delay_vector_ab,gain_vector_ab);

	%% 反射后，信号经过rayleigh信道
	rayleigh_chan_ab.ResetBeforeFiltering = 0;
	back_after_rayleigh_ab = filter(rayleigh_chan_ab,data_a_back); 

	%% Bob接收的信号为Alice反射信号和PT直传信号的叠加
	data_total_b = back_after_rayleigh_ab + data_after_rayleigh_b;

	%% 为叠加信号添加高斯噪声
	SNR_tag_b = 9;
	data_b = awgn(data_total_b,SNR_tag_b,'measured');

	%% 计算在Bob端接收信号的功率
	[Pxx_hamming_b, F_b] = periodogram(data_b,hamming(length(data_b)),[],Fs,'centered','psd');
	power_b = bandpower(Pxx_hamming_b,F_b,'psd');
	power_b_db = 10*log10(power_b/2);

	power_b_array(index) = power_b;

	%%**********************************************************************************************
	%% bob的反射路径，将信号反射给Alice
	coeffi_b = 0.9;						%% 反射因子
	data_b_coeffi = data_tag_b.*coeffi_b;

	%% 进行ASK调制，使用已经设置好的方波序列
	data_b_back = data_b_coeffi.*base_signal;

	%% 使用与Alice到Bob之间的信道模型
	% %% 构造rayleigh信道，表示Alice与Bob之间的信道
	% delay_vector_ab = [0, 4, 10, 19, 31, 43]*1e-9; 	% Discrete delays of four-path channel (s)
	% gain_vector_ab  = [0 -1.0 -5.0 -10.0 -18.0 -26.0]; 	% Average path gains (dB)			
	% max_Doppler_shift_ab = 10;  					% Maximum Doppler shift of diffuse components (Hz)			
	% rayleigh_chan_ab = rayleighchan(T,max_Doppler_shift_ab,delay_vector_ab,gain_vector_ab);

	%% 反射后，信号经过rayleigh信道
	rayleigh_chan_ab.ResetBeforeFiltering = 0;
	back_after_rayleigh_ba = filter(rayleigh_chan_ab,data_b_back); 

	%% Alice接收的信号为Bob反射信号和PT直传信号的叠加
	data_total_a = back_after_rayleigh_ba + data_after_rayleigh_a;

	%% 为叠加信号添加高斯噪声
	SNR_tag_a = 8;
	data_a = awgn(data_total_a,SNR_tag_a,'measured');


	%% 计算在Alice端接收信号的功率
	[Pxx_hamming_a, F_a] = periodogram(data_a,hamming(length(data_a)),[],Fs,'centered','psd');
	power_a = bandpower(Pxx_hamming_a,F_a,'psd');
	power_a_db = 10*log10(power_a/2);

	power_a_array(index) = power_a;

end

figure;
plot(num,power_tag_a_array);
hold on;
plot(num,power_tag_b_array);
hold on;
plot(num,power_b_array);
hold on;
plot(num,power_a_array);
legend('PT-a','PT-b','total-b','total-a');
grid on;


%% 计算序列的相关性
r = corr2(power_b_array,power_a_array)

alice = (power_a_array-power_tag_a_array).*power_tag_a_array;
bob = (power_b_array-power_tag_b_array).*power_tag_b_array;

cor = corr2(alice,bob)

figure;
plot(num,power_a_array.*power_tag_a_array);
hold on;
plot(num,power_b_array.*power_tag_b_array);
legend('alice','bob');
grid on;