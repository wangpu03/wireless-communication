%==========================================================================
%Name:      spectrum_analysis.m
%Desc:      以高斯脉冲信号为例，求解其频谱、双边功率谱、单边功率谱、双边功率谱密度、
%           单边功率谱密度，这里高斯信号的半波全宽FWHM=50ps，中心点位于2.5ns处。
%Parameter: 
%Return:    
%Author:    wangpu(wangpu03@gmail.com)
%Date:      2018-2-21
%高斯函数公式有两种表达方式 a*exp(-((x-b)/2c)^2)，和a*exp(-((x-b)/c)^2)
%其中c和2c都是反映脉冲陡度。对于前者，2ln2c表示波峰半高度时的宽度，后者则为2lnc

%fftshift是针对频域的，将FFT的DC分量移到频谱中心，操作就是对换数据的左右两边
%例子
%x=[1 2 3 4]
%fftshift(x) -> [3 4 1 2]
%直接用fft得出的数据与频率不是对应的，fftshift将其纠正过来
%如果期望绘制的幅频图的频率范围为0~fs，则无需运行fftshift变换
%如果期望绘制的幅频图的频率范围为-fs/2~fs/2，则需要运行fftshift变换
%=========================================================================

clc;
clear;
FWHM=50e-12;            %高斯信号FWHM宽度，为50ps
time_window=100*FWHM;   %高斯信号的采样窗口宽度，该值决定了傅里叶变换后的频率分辨率
Ns=2048;                %采样点
dt=time_window/(Ns-1);  %采样时间间隔
t=0:dt:time_window;     %采样时间

gauss_time=exp(-0.5*(2*sqrt(2*log(2))*(t-2.5e-9)/FWHM).^2); %高斯脉冲，中心位于2.5ns处。

plot(t*1e+9,gauss_time,'linewidth',0.5);
xlabel('Time/ns');
ylabel('Amplitude/V');
title('Gauss pulse');


%===========以下计算双边谱、双边功率谱、双边功率谱密度=================
gauss_spec=fftshift(fft(ifftshift(gauss_time)));    %傅里叶变换，并且进行fftshift移位操作。
gauss_spec=gauss_spec/Ns;   						%求实际的幅度值；
df=1/time_window;               					%频率分辨率
k=floor(-(Ns-1)/2:(Ns-1)/2);    
% k=0:Ns-1;
double_f=k*df;   									%双边频谱对应的频点


figure(2); %幅度谱
plot(double_f*1e-9,abs(gauss_spec),'linewidth',0.5);
xlabel('Frequency/GHz');
ylabel('Amplitude/V');
title('double Amplitude spectrum');


figure(3); %相位谱
plot(double_f*1e-9,angle(gauss_spec),'linewidth',0.5);
xlabel('Frequency/GHz');
ylabel('Phase/rad');
title('double Phase spectrum');


figure(4); %功率谱
double_power_spec_W=abs(gauss_spec).^2;                 %双边功率谱，单位W；
double_power_spec_mW=double_power_spec_W*1e+3;          %双边功率谱，单位mW；
double_power_spec_dBm=10*log10(double_power_spec_mW);   %双边功率谱，单位dBm；

plot(double_f*1e-9,double_power_spec_dBm,'linewidth',0.5);
xlabel('Frequency/GHz');
ylabel('Power/dBm');
title('double Power spectrum');


figure(5); %功率谱密度
double_power_specD_W=abs(gauss_spec).^2/(df);       	%双边功率谱密度,单位W/Hz
double_power_specD_mW=double_power_specD_W*1e+3;    	%双边功率谱密度,单位mW/Hz
double_power_specD_dBm=10*log10(double_power_specD_mW);	%双边功率谱密度,单位dBm/Hz

plot(double_f*1e-9,double_power_specD_dBm,'linewidth',0.5);
xlabel('Frequency/GHz');
ylabel('Power/(dBm/Hz)');
title('double power spectrum Density');


%==========以下计算单边谱、单边功率谱及单边功率谱密度=========
gauss_spec=fft(ifftshift(gauss_time));  %计算单边谱无需fftshift
gauss_spec=gauss_spec/Ns;       		%计算真实的幅度值
single_gauss_spec=gauss_spec(1:floor(Ns/2));
single_f=(0:floor(Ns/2)-1)*df;


figure(6); %幅度谱
plot(single_f*1e-9,abs(single_gauss_spec),'linewidth',0.5);
xlabel('Frequency/GHz');
ylabel('Amplitude/V');
title('single Amplitude spectrum');


figure(7); %相位谱
plot(single_f*1e-9,angle(single_gauss_spec),'linewidth',0.5);
xlabel('Frequency/GHz');
ylabel('Phase/rad');
title('single Phase spectrum');

figure(8);%功率谱
double_power_spec_W=abs(gauss_spec).^2;  
single_power_spec_W=2*double_power_spec_W(1:floor(Ns/2));   %单边功率谱，单位W
single_power_spec_mW=single_power_spec_W*1e+3;              %单边功率谱，单位mW；
single_power_spec_dBm=10*log10(single_power_spec_mW);       %双边功率谱，单位dBm；
plot(single_f*1e-9,single_power_spec_dBm,'linewidth',0.5);  
xlabel('Frequency/GHz');
ylabel('Power/dBm');
title('single Power spectrum');


figure(9);%功率谱密度
double_power_specD_W=abs(gauss_spec).^2/(df);
single_power_specD_W=2*double_power_specD_W(1:floor(Ns/2));         %单边功率谱密度，单位W/Hz
single_power_specD_mW=single_power_specD_W*1e+3;                    %单边功率谱密度，单位mW/Hz
single_power_specD_dBm=10*log10(single_power_specD_mW);             %单边功率谱密度，单位dBm/Hz

plot(single_f*1e-9,single_power_specD_mW,'linewidth',0.5);
xlabel('Frequency/GHz');
ylabel('Power/(dBm/Hz)');
title('single power spectrum density');