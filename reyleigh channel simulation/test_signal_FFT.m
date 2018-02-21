%% x=0.5*sin(2*pi*15*t)+2*sin(2*pi*40*t), 采样频率fs=100Hz，分别绘制N=128、1024点幅频图
% fs=100Hz，Nyquist频率为fs/2=50Hz。整个频谱图是以Nyquist频率为对称轴的。
% 并且可以明显识别出信号中含有两种频率成分：15Hz和40Hz。
% 由此可以知道FFT变换数据的对称性。因此用FFT对信号做谱分析，
% 只需考察0~Nyquist频率范围内的福频特性。若没有给出采样频率和采样间隔，
% 则分析通常对归一化频率0~1进行。另外，振幅的大小与所用采样点数有关，
% 采用128点和1024点的相同频率的振幅是有不同的表现值，但在同一幅图中，
% 40Hz与15Hz振动幅值之比均为4：1，与真实振幅0.5：2是一致的。
% 为了与真实振幅对应，需要将变换后结果乘以2除以N。


%对信号采样数据为128点的处理
clf;
fs = 100;
N = 128;			%采样频率和数据点数
n = 0:N-1;
t = n/fs;			%时间序列

x=0.5*sin(2*pi*15*t)+2*sin(2*pi*40*t); %信号

y = fft(x,N);		%对信号进行快速fourier变换
mag = abs(y);		%求得Fourier变换后的幅值

f = n*fs/N;			%频率序列

subplot(2,2,1);		%绘出随频率变化的振幅
plot(f,mag*2/N);
xlabel('Frequency/Hz');
ylabel('Amplitude');title('N=128');grid on;

subplot(2,2,2);
plot(f(1:N/2),mag(1:N/2)*2/N); %绘出Nyquist频率之前随频率变化的振幅
xlabel('Frequency/Hz');
ylabel('Amplitude');title('N=128');grid on;


%对信号采样数据为1024点的处理
fs=100;
N=1024;
n=0:N-1;
t=n/fs;

x=0.5*sin(2*pi*15*t)+2*sin(2*pi*40*t); %信号

y=fft(x,N);    %对信号进行快速Fourier变换
mag=abs(y);    %求取Fourier变换的振幅
f=n*fs/N;

subplot(2,2,3),plot(f,mag*2/N); %绘出随频率变化的振幅
xlabel('Frequency/Hz');
ylabel('Amplitude');title('N=1024');grid on;

subplot(2,2,4)
plot(f(1:N/2),mag(1:N/2)*2/N); %绘出Nyquist频率之前随频率变化的振幅
xlabel('Frequency/Hz');
ylabel('Amplitude');title('N=1024');grid on;