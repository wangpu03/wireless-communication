%% 相位谱
% 例子： x = sin(2*pi*60*t);
clc;
fs = 1024;
N = 10240;
n = 0:N-1;
t = n/fs;
x = sin(2*pi*60*t);

y = fft(x,N);
A = abs(y);
f = n*fs/N;
ph = 2*angle(y(1:N/2));

ph = ph*180/pi;
subplot(2,1,1)
plot(f(1:N/2),ph(1:N/2));
xlabel('F/Hz');
ylabel('phase');
title('phase spectrum');
grid on;

