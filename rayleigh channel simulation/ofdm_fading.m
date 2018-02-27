%   system.m
%   主程序

%% 仿真信号在室内，经过大尺度和小尺度衰落后的OFDM调制

clear all
clc
t0=clock;
format long
cp_length=16;%cp长度为16
N_carrier=64;%802.11a OFDM子载波个数64
bandwidth=17000000;%系统带宽为17MHz
t_interval=(1/bandwidth)*N_carrier/(cp_length+N_carrier);%采样间隔64/70us，加上循环前缀后，采样率增加
delta_f=bandwidth/N_carrier;%Hz
sendpower=800;%mW
d= 10;                         %传输距离m
w= [0]  ;           %穿过的墙壁损耗列向量db
f= [0] ;                     %穿过的地板损耗列向量db
SNR_dB=[0 4 8 12 16 20 24 28 32];%Eb/N0
SJR_dB=[0];
row_num=length(SNR_dB);
ber_snr_persjr_qpsk=zeros(row_num,length(SJR_dB));%最终给出各种SJR下误码率随SNR变化的曲线
ber_snr_persjr_16qam=zeros(row_num,length(SJR_dB));
counter_d=1000000;%各loop信道抽样跳动间隔
no=[25,25,25,25,25,25];
for jj=1:length(SJR_dB)%每个SJR得到一条曲线
    sjr=10^(SJR_dB(jj)/10);
    for i=1:length(SNR_dB)%每个SNR点上仿真若干次
    error_bit=0;%错误比特数统计
    total_bit_num=0;%发送总比特数统计
    error_bit16=0;%错误比特数统计
    total_bit_num16=0;%发送总比特数统计
    loop_num=5; %共仿真5次
    for l=1:loop_num
    ofdm_symbol_num=12;%每次仿真产生200个ofdm符号,则每次仿真共有200×64个星座映射符号；QPSK调制下，1个星座映射符号包含2个bit
    map_flag=4;
    code_rate=1/2;
    bit_source= randn(1,N_carrier*ofdm_symbol_num*map_flag) > 0;
    [nbit,mbit]=size(bit_source);
    total_bit_num=total_bit_num+nbit*mbit;
    %以下为发送端卷积码编码和交织操作
    % Convolutional encoding
    coded_bit_stream = tx_conv_encoder(bit_source);
    tx_bits = tx_puncture(coded_bit_stream,code_rate);
    rdy_to_mod_bits = tx_make_int_num_ofdm_syms(tx_bits,N_carrier,map_flag);
    rdy_to_mod_bits1 =tx_interleaver(rdy_to_mod_bits,N_carrier,map_flag);
    rdy_to_mod_bits2 = reshape(rdy_to_mod_bits1,N_carrier*map_flag,length(rdy_to_mod_bits1)/(N_carrier*map_flag));%串行比特流转变成比特矩阵
      map_out=map_module(rdy_to_mod_bits2,map_flag); 
      ofdm_modulation_out=sqrt(N_carrier)*ifft(map_out,N_carrier);%作64点逆FFT运算，完成ofdm调制,前面乘系数sqtr(64)是为了保持ifft前后的符号能量不变
        ofdm_cp_out=insert_cp(ofdm_modulation_out,cp_length);%插入循环前缀 
        snr=10^(SNR_dB(i)/10);
    [nnl,mml]=size(ofdm_cp_out);

 %*******************大尺度衰落***********************%
    receivepower=bigfade(sendpower,d,w,f);
    sgma=sqrt(receivepower*t_interval/(2*snr)/map_flag);

 %************** 以下过程为ofdm符号通过频率选择性多径信道********% 
    num=6;
    %假设功率延迟谱服从负指数分布~exp(-t/trms),trms=(1/4)*cp时长；
    %t在0~cp时长上均匀分布
    %若cp时长为16e-6s，可以取5径延迟如下
    delay=[0 5e-8 11e-8 17e-8 29e-8 31e-8 ];
    trms=1.5e-6;
    var_pow=[0 -3 -10 -18 -26 -32];%各径功率衰减,以dB形式给出
    fd=50;%最大doppler频率为50Hz
     counter_begin=(l-1+1000000)*5*counter_d;
    %信道采样点数，每个调制符号采一个点
    ofdm_cp_out1=ofdm_cp_out*sqrt(receivepower*t_interval);
     [passchan_ofdm_symbol,Hk]=multipath_chann(ofdm_cp_out1,num,var_pow,delay,fd,t_interval,counter_begin,counter_d,cp_length,no);%小尺度
     cutcp_ofdm_symbol1=cut_cp(passchan_ofdm_symbol,cp_length);
     ofdm_demodulation_out1=fft(cutcp_ofdm_symbol1,N_carrier)/sqrt(N_carrier);
     HHk=ofdm_demodulation_out1./map_out;
    %******** 以上过程为ofdm符号通过频率选择性多径信道********

    %******** 以下过程为ofdm符号加高斯白噪声 **************
    passnoise_ofdm_symbol=add_noise(sgma,passchan_ofdm_symbol);%加入随机高斯白噪声，receive_ofdm_symbol为最终接收机收到的ofdm符号块
    
    cutcp_ofdm_symbol=cut_cp(passnoise_ofdm_symbol,cp_length);%去除循环前缀
    ofdm_demodulation_out=fft(cutcp_ofdm_symbol,N_carrier)/sqrt(N_carrier);%作128点FFT运算，完成ofdm解调
    ofdm_demodulation_out=ofdm_demodulation_out./HHk;
    receive_sig=de_map_module(ofdm_demodulation_out,map_flag);%解映射
     %以下为接收端解交织和卷积码译码操作
     deint_bits = rx_deinterleave(receive_sig,N_carrier,map_flag);
     depunc_bits = rx_depuncture(deint_bits,code_rate);
     viterbi_input = depunc_bits(1:(N_carrier*ofdm_symbol_num*map_flag+6)*2);
     % Vitervi decoding
     data_bits = rx_viterbi_decode(viterbi_input);
     receive_bit_sig=data_bits(1:N_carrier*ofdm_symbol_num*map_flag);
    %以下过程统计接收信号中的错误比特数
    err_num=error_count(bit_source,receive_bit_sig);
    error_bit=error_bit+err_num;
    fprintf('error_bit is %d\n',error_bit);
 end
 % %计算各种信噪比下的误比特率
ber_snr_persjr_qpsk(i,jj)=error_bit/total_bit_num;
fprintf('ber_snr_persjr_qpsk is %f\n',ber_snr_persjr_qpsk(i,jj));
end
end
   elapse_time=etime(clock,t0);
   semilogy(SNR_dB,ber_snr_persjr_qpsk);


% bigfade.m
% 大尺度衰落
% Generate bigscale fading
function[receivepower]=bigfade(sendpower,d,w,f)

%****************** variables *********************
L= 61.6 ;           %第一米路径损耗db,2.4G时为40
n=4;                 %路径损耗指数
%**************************************************

Lptco=L+10*n*log10(d)+sum(w)+sum(f);
s=10^(-Lptco/10);
receivepower=s*sendpower;
%******************** end of file *****************



% multipath_chann.m
function [output_sig,Hk]=multipath_chann(input_sig,num,var_pow,delay,fd,t_interval,counter_begin,counter,cp_n,no)
%input_sig输入信号矩阵,加了cp后的信号，大小为NL×(子载波个数+cp长度lp)；
%num多径数;
%var_pow各径相对主径的平均功率,单位dB；
%delay各径延时,单位s；
%fd最大dopple频率；
%t_interval为离散信道抽样时间间隔，等于OFDM符号长度/(子载波个数+cp长度lp)；
%counter各径间隔记录
%count_begin本次产生信道开始记录的初始位置
%no各径用于叠加形成Raleigh衰落的正弦波条数
t_shift=round(delay/t_interval);%归一化各径延时
[nl,l]=size(input_sig);
output_sig=zeros(size(input_sig));
chann_l=nl*l;%信道采样点数，若一个调制符号采样一个信道点，则信道采样点数等于输入信号中的调制符号个数
selec_ray_chan=zeros(num,chann_l);%初始化频率选择性信道，径数＝num
pow_per_channel=10.^(var_pow/10);%各径功率线性化，从dB转变成线性
total_pow_allchan=sum(pow_per_channel(1:num));%各径功率之和
%以下for循环产生相互独立的num条rayleigh信道
for k=1:num
    atts=sqrt(pow_per_channel(k));
    selec_ray_chan(k,:)=atts*rayleigh_fade(chann_l,t_interval,fd,counter_begin+k*1000,no(k))/sqrt(total_pow_allchan);
  end
%以下计算信道频率响应值，行数与载波数相同，列数与符号个数相同
N_carrier=nl-cp_n;
N_ofdm=l;
selec_ray_channew=zeros(num,chann_l+N_ofdm*(N_carrier+cp_n));
selec_ray_channew(:,1:chann_l)=selec_ray_chan;
Hk=zeros(N_carrier,N_ofdm);
for kk=1:N_ofdm
    for ii=1:N_carrier
        sum1=0;
        for ll=1:num
            sum1=sum1+selec_ray_channew(ll,(kk-1)*(N_carrier+cp_n)+cp_n+ii+t_shift(ll))*exp(-sqrt(-1)*2*pi*(ii-1)*t_shift(ll)/N_carrier);
        end
        Hk(ii,kk)=sum1;
    end
end
for k=1:l
    input_sig_serial(((k-1)*nl+1):k*nl)=input_sig(:,k).';%输入信号矩阵转变成串行序列
end
delay_sig=zeros(num,chann_l);%初始化延时后的送入各径的信号，每径所含符号数为chann_l
%以下for循环为各径的输入信号做延迟处理
for f=1:num
    if t_shift(f)~=0
        delay_sig(f,1:t_shift(f))=zeros(1,t_shift(f));
    end
    delay_sig(f,(t_shift(f)+1):chann_l)= input_sig_serial(1:(chann_l-t_shift(f)));
end
output_sig_serial=zeros(1,chann_l);%初始化输出信号串行序列
%得到各径叠加后的输出信号序列
for f=1:num
        output_sig_serial= output_sig_serial+selec_ray_chan(f,:).*delay_sig(f,:);
end
for k=1:l
    output_sig(:,k)=output_sig_serial(((k-1)*nl+1):k*nl).';%输出信号串行序列转变成与输入信号相同的矩阵形式，做为本函数输出
end


% rayleigh_fade.m
function ray_chann=rayleigh_fade(nsamp,tstp,fd,counter,no)
%****************** variables *************************
% idata  : input Ich data     
% qdata  : input Qch data     
% iout   : output Ich data
% qout   : output Qch data
% ramp   : Amplitude contaminated by fading
% rcos   : Cosine value contaminated by fading
% rsin   : Cosine value contaminated by fading
% nsamp  : Number of samples to be simulated       
% tstp   : Minimum time resolution                    
% fd     : maximum doppler frequency               
% no     : number of waves in order to generate fading   
% counter  : fading counter                          
% flat     : flat fading or not 
% (1->flat (only amplitude is fluctuated),0->nomal(phase and amplitude are fluctutated)    
%******************************************************
% no=25;
if fd ~= 0.0  
    ac0 = sqrt(1.0 ./ (2.0.*(no + 1)));   % power normalized constant(ich)
    as0 = sqrt(1.0 ./ (2.0.*no));         % power normalized constant(qch)
    pai = 3.14159265;   
    wm = 2.0.*pai.*fd;
    n = 4.*no + 2;
    ts = tstp;
    wmts = wm.*ts;
    paino = pai./no;                        

    xc=zeros(1,nsamp);
    xs=zeros(1,nsamp);
    ic=[1:nsamp]+counter;

  for nn = 1: no
	  cwn = cos( cos(2.0.*pai.*nn./n).*ic.*wmts );
	  xc = xc + cos(paino.*nn).*cwn;
	  xs = xs + sin(paino.*nn).*cwn;
  end
  cwmt = sqrt(2.0).*cos(ic.*wmts);
  xc = (2.0.*xc + cwmt).*ac0;
  xs = 2.0.*xs.*as0;
  ray_chann=sqrt(xc.^2+xs.^2);
else  
 ray_chann=ones(1,nsamp);
end

% ************************end of file***********************************
