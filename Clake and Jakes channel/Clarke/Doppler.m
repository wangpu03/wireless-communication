function y=Doppler(fd,Nfft)
% fd 最大多普勒频率
% Nfft 频域点数
df=2*fd/Nfft;
%第一个DC为0
f(1)=0;
y(1)=1.5/(pi*fd);
%单侧普的其他分量
for i=2:Nfft/2
    f(i)=(i-1)*df;
    y([i Nfft-i+2])=1.5/(pi*fd*sqrt(1-(f(i)/fd)^2));
end
% s使用最后3个采样点实现多项式拟合
nFitPoints=3;
kk=[Nfft/2-nFitPoints:Nfft/2];
polyFreq=polyfit(f(kk),y(kk),nFitPoints);
y(Nfft/2+1)=polyval(polyFreq,f(Nfft/2)+df);
