function [h,tf]=Jakes_model(fd,Ts,Ns,t0,E0,phi_N)
% inputs:
%       fd:  doppler frequency
%       Ts: 采样周期
%       Ns: 采样点数
%       t0: 初始时间
%       E0: 信道功率
%       phi_N: 具有最大多普勒频率信号的初始相位
% output:
%       h: 复合衰落信道
%       t_state: 当前时刻

% check inputs
if nargin<6, phi_N=0; end
if nargin<5, E0=1; end
if nargin<4, t0=0; end
if nargin<3
    error('More inputs are needed.');
end

N0=8;   %
N=4*N0+2;   % 多径数
wd=2*pi*fd;
t=t0+[0:Ns-1]*Ts;
tf=t(end)+Ts;   % the end time
% 尽量利用矩阵
% hI=2*sum(cos(phi_n)*coswt)+sqrt(2)*cos(phi_N)*cos(wd*t)
% hQ=2*sum(sin(phi_n)*coswt)+sqrt(2)*sin(phi_N)*cos(wd*t)
coswt=[sqrt(2)*cos(wd*t).*cos(phi_N);2*cos(wd*cos(2*pi/N*[1:N0]')*t)];
h=E0/sqrt(2*N0+1)*exp(j*[phi_N pi/(N0+1)*[1:N0]])*coswt;

