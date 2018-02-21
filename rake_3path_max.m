%Rake reciever/rake.m
%
clear all;
close all;
Tc = 1;
N = 128;
gx = '1000010111000101';%g(x) = x^15+x^13+x^9+x^8+x^7+x^5+1
g = bin2dec(gx);
state = 1;
L = 2^13;

EcN0dB = -21:-14;   %different noise value

for  k = 1:length(EcN0dB)
    error(k) = 0;%the number of error bits
    total(k) = 0;%the total number of transmission bits
    sigma(k) = sqrt(10.^(-EcN0dB(k)/10)/2);
    while(error(k)<100)
        %3 channels
        p1 = sqrt(0.5/2)*(randn(1,L)+j*randn(1,L)); 
        p2 = sqrt(0.3/2)*(randn(1,L)+j*randn(1,L));
        p3 = sqrt(0.2/2)*(randn(1,L)+j*randn(1,L));
        %corresponding relay times
        t1 = 0;
        t2 = 1;
        t3 = 2;
        
        pt = mgen(g,state,L+t3);%call the M sequence generation function
        pt = 2*pt-1;
        
        %generate 64 random bits (-1,0,1)
        d = sign(randn(1,L/N));
        
        %insert the 127 0s in the sequence and expand 64 bits to 128*64
        dd = sigexpand(d,N);                
        s = conv(dd,ones(1,N));             %convolution
        st = s(1:L+t3).*pt(1:L+t3);         %spread frequency
        
        %add the noise
        z = sigma(k)*(randn(1,L)+j*randn(1,L)); 
        
        rt = st(1:L).*p1 + st(t2+1:L+t2).*p2 + st(t3+1:t3+L).*p3 + z;
        
        
        %rake receiving
        r1 = rt.*conj(p1).*pt(1:L);
        r2 = rt.*conj(p2).*pt(t2+1:L+t2);
        r3 = rt.*conj(p3).*pt(t3+1:L+t3);
        
        %integration
        r1 = reshape(r1,N,L/N);y1 = sum(r1);
        r2 = reshape(r2,N,L/N);y2 = sum(r2);
        r3 = reshape(r3,N,L/N);y3 = sum(r3);
        
        %sum
        y = y1 + y2 + y3;%?????
        
        %decode
        dc = sign(real(y));
        
        error(k) = error(k) + sum(abs((d-dc))/2);
        total(k) = total(k) + L/N;
        BitErrorRate(k) = error(k)/total(k);
    end
end

%draw  the curved line about the bit error rate
semilogy(EcN0dB,BitErrorRate);
grid on;