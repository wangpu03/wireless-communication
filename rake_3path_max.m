%Rake接收机rake.m 
clear all; 
close all; 
Tc = 1; 
N = 128; 
gx = '1000010111000101';%g(x) = x^15+x^13+x^9+x^8+x^7+x^5+1 
g = bin2dec(gx); 
state = 1; 
L = 2^13; 
 
EcN0dB = -21:-14; 
 
for  k = 1:length(EcN0dB) 
    error(k) = 0;%计数错误比特数 
    total(k) = 0;%计数总的传输比特数 
    sigma(k) = sqrt(10.^(-EcN0dB(k)/10)/2); 
    while(error(k)<100) 
        %多径结构 
        p1 = sqrt(0.5/2)*(randn(1,L)+j*randn(1,L)); 
        p2 = sqrt(0.3/2)*(randn(1,L)+j*randn(1,L)); 
        p3 = sqrt(0.2/2)*(randn(1,L)+j*randn(1,L)); 
        t1 = 0; 
        t2 = 1; 
        t3 = 2; 
         
        pt = mgen(g,state,L+t3);%调用m序列发生器函数 
        pt = 2*pt-1; 
         
        %数据产生 
        d = sign(randn(1,L/N));%一次64个 
        %扩频，先将数据扩展，然后与pt点积 
        dd = sigexpand(d,N); 
        s = conv(dd,ones(1,N)); 
        st = s(1:L+t3).*pt(1:L+t3);         %扩频 
         
        %经过多径信道，加入噪声 
        z = sigma(k)*(randn(1,L)+j*randn(1,L)); 
         
        rt = st(1:L).*p1 + st(t2+1:L+t2).*p2 + st(t3+1:t3+L).*p3 + z; 
        %rake接收 
        r1 = rt.*conj(p1).*pt(1:L); 
        r2 = rt.*conj(p2).*pt(t2+1:L+t2); 
        r3 = rt.*conj(p3).*pt(t3+1:L+t3); 
         
        %积分 
        r1 = reshape(r1,N,L/N);y1 = sum(r1); 
        r2 = reshape(r2,N,L/N);y2 = sum(r2); 
        r3 = reshape(r3,N,L/N);y3 = sum(r3); 
        %合并 
        y = y1 + y2 + y3;%最大比合并 
        %判决 
        dc = sign(real(y)); 
        error(k) = error(k) + sum(abs((d-dc))/2); 
        total(k) = total(k) + L/N; 
        BitErrorRate(k) = error(k)/total(k); 
    end 
end 
semilogy(EcN0dB,BitErrorRate); 
grid on;

