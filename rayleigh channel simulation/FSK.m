clear  
clc  
  
%% 参数设置  
Modulation_Order=2;  
Packet_Length=16;  
SymbolRate=1000;  
SamplesPerSymbol=8;%每个符号的样点数  
df=400;%FSK调制相邻载波频率间隔，单位Hz  
fs=SamplesPerSymbol*SymbolRate;  
SNR=10;  
ts=1/fs;  
fd=0;  
k=3;  
tau=[0 1 2]*ts;%每条径的时延  
pdb=[2 1 0];%每条径的增益  
  
  
%% 调制器、解调器、莱斯信道设置  
H_Mod = comm.FSKModulator('ModulationOrder',Modulation_Order,'SymbolRate',1000 ,'SamplesPerSymbol',8, 'FrequencySeparation',df); %默认连续相位、整数输入   
H_Demod = comm.FSKDemodulator('ModulationOrder',Modulation_Order,'SymbolRate',1000 ,'SamplesPerSymbol',8, 'FrequencySeparation',df);  
chan_Rici = ricianchan(ts,fd,k,tau,pdb);  
  
%% 调制 过信道 加噪声 解调  
data_Original = randi([0 Modulation_Order-1],Packet_Length,1);  
data_AfterFSK=step(H_Mod,data_Original);  
data_AfterRici=filter(chan_Rici,data_AfterFSK);  
data_AfterAwgn=awgn(data_AfterRici,SNR,'measured');  
data_AfterDeFSK=step(H_Demod,data_AfterAwgn);  
  
%% 画图  
  
% FSK调制后的信号为实际信号的等效低通信号，即复数信号，看波形时只要看实部即可。  
figure(1)  
plot(real(data_AfterFSK),'r*-');  
  
figure(2)  
plot(data_Original,'r*-');  
hold on  
plot(data_AfterDeFSK,'go-');  