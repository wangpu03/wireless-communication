%% Theoretical Power of Multiple Sinusoids

clc;
Fs = 1024;
t = 0:1/Fs:1-(1/Fs);

A0 = 1.5;					% Vpeak, DC
A1 = 4;
A2 =3;

F1 = 100;
F2 = 200;

x = A0 + A1*sin(2*pi*F1*t) + A2*sin(2*pi*F2*t);

% let's look at a portion of our signal
idx = 1:200;
plot(t(idx),x(idx));
grid on;
ylabel('Amplitude');
xlabel('Time(sec)');

hgcf = gcf;
hgcf.Color = [1 1 1];

%% The theoretical average power (mean-square) of each complex sinusoid is A^2/4
%% the theoretical average power of each complex sinusoid
power_theoretical = A0^2 + (A1^2/4)*2 + (A2^2/4)*2

% use the periodogram function to calculate and plot the power spectrum of the signal.
periodogram(x, hamming(length(x)),[],Fs,'centered','power');

% estimating the signal's total average power by "integrating" under the PSD curve we get:
% 方括号中，前者是功率，后者是频率，构成功率谱
[Pxx, F] = periodogram(x, hamming(length(x)), [], Fs, 'centered', 'psd');

power_freqdomain = bandpower(Pxx,F,'psd');

power_timedomain = sum(abs(x).^2)/length(x);

%% 