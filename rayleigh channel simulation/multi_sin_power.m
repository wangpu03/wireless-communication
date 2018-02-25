%% Theoretical Power of Multiple Sinusoids

clc;
clear;
close all;
Fs = 1024;
t = 0:1/Fs:1-(1/Fs);

A0 = 1.5;					% Vpeak, DC
A1 = 4;
A2 =3;

F1 = 100;
F2 = 200;

x = A0 + A1*sin(2*pi*F1*t) + A2*sin(2*pi*F2*t);

% let's look at a portion of our signal
figure(1);
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

figure(2);
% use the periodogram function to calculate and plot the power spectrum of the signal.
periodogram(x, hamming(length(x)),[],Fs,'centered','power');

% estimating the signal's total average power by "integrating" under the PSD curve we get:
% 方括号中，前者是功率，后者是频率，构成功率谱
[Pxx, F] = periodogram(x, hamming(length(x)), [], Fs, 'centered', 'psd');

power_freqdomain = bandpower(Pxx,F,'psd')

power_timedomain = sum(abs(x).^2)/length(x)

%% Relationship between Power Spectrum, Power Spectral Density and ENBW
Pxx = periodogram(x, hamming(length(x)), [], Fs, 'centered', 'psd');
Sxx = periodogram(x, hamming(length(x)), [], Fs, 'centered', 'power');

plot(F, Sxx ./ Pxx)
grid on
axis tight
xlabel('Frequency (Hz)')
title('Ratio between Power Spectrum and Power Spectral Density')

ratio = mean(Sxx ./ Pxx)

bw = enbw(hamming(length(x)),Fs);

%% Enhanced Power Measurements Using Reassigned Periodogram
%% In the previous sections, power was measured from one or multiple 
%% sinusoids having a frequency that coincided with a bin. Peak power 
%% estimates are usually less accurate when the signal frequency is out 
%% of bin. To see this effect, create a sinusoid with a non-integer number 
%% of cycles over a one second period.


Fs = 1024;
t  = 0:1/Fs:1-(1/Fs);
A = 1;
F = 20.4;
x = A*sin(2*pi*F*t);
NFFT = length(x);
power_theoretical = 10*log10((A^2/4)*2);

% Create a Hamming window and a flat top window.
w1 = hamming(length(x));
w2 = flattopwin(length(x));

% Compute the periodogram of x using the Hamming window. Zoom in on the peak.
h1 = figure; 
hold on
stem(F,power_theoretical,'BaseValue',-50);
[Pxx1,f1] = periodogram(x, w1, NFFT, Fs, 'power');
plot(f1,10*log10(Pxx1))
axis([0 40 -45 0])
legend('Theoretical','Periodogram')
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
title('Periodogram Power Spectrum Estimate')
grid on

% 取最大值
[Pmax,imax] = max(Pxx1);
dPmax_w1 = 10*log10(Pmax) - power_theoretical

% 最大值与中心频率的偏差
dFreq = f1(imax) - F

%% Reduce Amplitude Error with Zero-Padding
%% To see why this is happening, compute the periodogram using a larger number of FFT bins.
figure
hold on
periodogram(x, w1, 100*NFFT, Fs, 'power')
ax = gca;
ax.ColorOrderIndex = 2;
plot(f1,10*log10(Pxx1),'+')
stem(F,power_theoretical,'BaseValue',-50)
axis([0 40 -40 0])
legend('NFFT = 1024','NFFT = 102400','Theoretical Peak')


figure(h1)
[Pxx,F1] = periodogram(x, w2, NFFT, Fs, 'power');
plot(F1,10*log10(Pxx))
legend('Theoretical','Hamming','Flat Top')

% 对比三种数值，其中两个已经重合，一模一样
[RPxx1,~,~,Fc1] = periodogram(x, w1, NFFT, Fs, 'power','reassigned');
[RPxx2,~,~,Fc2] = periodogram(x, w2, NFFT, Fs, 'power','reassigned');
figure
hold on
stem(F,power_theoretical,'BaseValue',-40)
stem(Fc1,10*log10(RPxx1),'BaseValue',-50)
stem(Fc2,10*log10(RPxx2),'BaseValue',-50)
legend('Theoretical','Hamming Reassignment','Flattop Reassignment')
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
title('Periodogram Power Spectrum Estimate')
axis([19.5 21 -4 -2])
grid on


[RPxx1max,imax1] = max(RPxx1);
[RPxx2max,imax2] = max(RPxx2);
dPmax_reassign_w1 = 10*log10(RPxx1max) - power_theoretical

dPmax_reassign_w2 = 10*log10(RPxx2max) - power_theoretical

Fc1(imax1)-F

Fc2(imax2)-F