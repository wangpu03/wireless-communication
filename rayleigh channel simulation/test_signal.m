%%  构造离散时间向量
Fs = 1000; 			% sampling frequency, 采用频率
T = 1/Fs;  			% sampling period, 采用周期
L = 1000; 			% length of signal
t = (0:L-1)*T; 		% time vector


f1 = 50;
f2 = 120;
%% 构造信号
%  信号包括一个振幅为0.7，频率为50Hz的正弦信号，一个振幅为1，频率为120的正弦信号
S = 0.7*sin(2*pi*f1*t)+sin(2*pi*f2*t);

%% 加噪声
X = S + 2*randn(size(t));


%% spike signal
N = 100;
s = zeros(N,1);
k = [20, 45, 70];
a = [2, -1, 1];
s(k) = a;

%% 四点脉冲信号
L1 = 4;
h = ones(L1,1)/4