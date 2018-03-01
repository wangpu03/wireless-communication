% num = 100;
% index = 1:num;
% for i = index
% 	disp(i);
% end


f1 = 900e6;				% 信号频率900MHz,信号周期为10/9 ns
f2 = 950e6;				% 信号频率900MHz,信号周期为10/9 ns
N  = 20;				% 信号周期内的采样点数
Fs = N*f1;				% sampling frequency, 采样频率
T = 1/Fs;  				% sampling period, 采样周期
L = 20*N; 			% length of signal

t = (0:L-1)*T;			% 采样时间s，fs的值越大，出来的波形失真越小
A = 4;					% 信号幅值


phase1 = 2*pi*rand();

%% 构造初始信号
source1 = A*sin(2*pi*f1*t);
[Pxx_hamming_b, F_b] = periodogram(source1,hamming(length(source1)),[],Fs,'centered','psd');
power_b = bandpower(Pxx_hamming_b,F_b,'psd');
power_b_db = 10*log10(power_b/2)

phase2 = 2*pi*rand();

%% 构造初始信号
source2 = A*sin(2*pi*f1*t);

y = source1 + source2;

[Pxx_hamming_a, F_a] = periodogram(y ,hamming(length(y )),[],Fs,'centered','psd');
power_a = bandpower(Pxx_hamming_a,F_a,'psd');
power_a_db = 10*log10(power_a/2)


figure;
plot(t,source1);
hold on;
plot(t,source2);

grid on;

plot(t,y);