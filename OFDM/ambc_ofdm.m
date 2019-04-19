%% simulation of the signals received by BDs over OFDM carrier
% but we do not simulate the backscattered model. only calculate the
% direct signal and backscatter signal without backsctter modulation

clear;
clc;

%% simulation parameters
% modulation methods: BPSK, QPSK,16QAM, 32QAM,64QAM
mod_method = 'QPSK';

% IFFT/FFT size
n_fft = 128;

% size of cyclic prefix extension
n_cpe = 32;

% target SNR (dB)
snr =20;

% number of channel taps (1 = no channel) 
n_taps1 = 8;        %f1
n_taps2 = 6;        %f1
n_taps_h12 = 3;        %h12

% channel estimaition method: none, LS
ch_est_method = 'LS';

% option to save plot to file
save_file = 0;

% calculate modulation order from modulation method
mod_methods = {'BPSK', 'QPSK','8PSK','16QAM', '32QAM','64QAM'};
mod_order = find(ismember(mod_methods, mod_method));

%% input data to binary stream
rand_ints = randi(2,25600,1)-1;
rand_bits = dec2bin(rand_ints(:));
im_bin = rand_bits(:);
%% binary stream to symbols 
% Parse binary stream into mod_order bit symbols
% pads input signal to appropriate length
sym_rem = mod(mod_order-mod(length(im_bin),mod_order),mod_order);
padding = repmat('0',sym_rem,1);
im_bin_padded = [im_bin;padding];
cons_data = reshape(im_bin_padded, mod_order,length(im_bin_padded)/mod_order)';
cons_sym_id = bin2dec(cons_data);

%% symbol modulation
% BPSK
if mod_order == 1
    mod_ind = 2^(mod_order-1);
    n = 0:pi/mod_ind:2*pi-pi/mod_ind;
    in_phase = cos(n);
    quadrature = sin(n);
    symbol_book = (in_phase+quadrature*1i)';
end

% phase shift keying about unit circle
if mod_order == 2 || mod_order ==3
    mod_ind = 2^(mod_order-1);
    n = 0:pi/mod_ind:2*pi-pi/mod_ind;
    in_phase = cos(n+pi/4);
    quadrature = sin(n+pi/4);
    symbol_book = (in_phase+quadrature*1i)';
end

% 16QAM, 64QAM modulaiton
if mod_order ==4 || mod_order ==6
     mod_ind = 2^(mod_order-2);
     in_phase = repmat(linspace(-1,1,mod_ind),mod_ind,1);
     quadrature =  repmat(linspace(-1,1,mod_ind)',1,mod_ind);
     symbol_book = in_phase(:)+quadrature(:)*1i;
end

% 32 QAM modulation
% generate 6x6 constellation and removes corners
if mod_order ==5
    mod_ind = 6;
    in_phase = repmat(linspace(-1,1,mod_ind),mod_ind,1);
    quadrature =  repmat(linspace(-1,1,mod_ind)',1,mod_ind);
    symbol_book = in_phase(:)+quadrature(:)*1i;
    symbol_book = symbol_book([2:5 7:30 32:35]);
end

% modulate data according to symbol_book
X = symbol_book(cons_sym_id+1);

%% use IFFT to move to time domain
% pad input signal to appropriate length
fft_rem = mod(n_fft-mod(length(X),n_fft),n_fft);
X_padded = [X;zeros(fft_rem,1)];
X_blocks = reshape(X_padded, n_fft,length(X_padded)/n_fft);
x = ifft(X_blocks);

% add cyclic prefix extension and shift from parallel to serial
x_cpe = [x(end-n_cpe+1:end,:);x];
x_s = x_cpe(:);

%*********************************************************************

%% add awgn
% calculation data power
data_pwr = mean(abs(x_s).^2);

% calculate the cp data power for n frame OFDM
x_s_cp = x_s(9:32);
for n = 1:1:9
    x_s_cp =[x_s_cp(:)' x_s(160*n+9:160*n+32)']';
end
data_cpe_pwr = mean(abs(x_s_cp).^2);

% add noise to channel 
noise_pwr = data_pwr/10^(snr/10);
noise = normrnd(0,sqrt(noise_pwr/2),size(x_s))+normrnd(0,sqrt(noise_pwr/2),size(x_s))*1i;
%x_s_noise = x_s+noise;

x_s_noise = x_s;

% measure SNR
snr_meas = 10*log10(mean(abs(x_s).^2)/mean(abs(noise).^2));

%% apply fading channel
% there is a big difference with or without the parameter 'same'
% with the 'same' returns only the central part of the convolution

%*********************************************************************
% f1
g1 = exp(-(0:n_taps1-1));
g1 = g1/norm(g1);
%x_s_noise_fading = conv(x_s_noise,g,'same');

x_s_noise_fading1 = conv(x_s_noise,g1);
%*********************************************************************
% h21
g21 = exp(-(0:n_taps_h12-1));
g21 = g21/norm(g21);
%x_s_noise_fading = conv(x_s_noise,g,'same');

x_s_noise_fading1_h21 = conv(x_s_noise_fading1,g21);
%*********************************************************************
% f2
g2 = exp(-(0:n_taps2-1));
g2 = g2/norm(g2);
%x_s_noise_fading = conv(x_s_noise,g,'same');

x_s_noise_fading2 = conv(x_s_noise,g2);
%*********************************************************************
% h12
g12 = exp(-(0:n_taps_h12-1));
g12 = g12/norm(g12);
%x_s_noise_fading = conv(x_s_noise,g,'same');

x_s_noise_fading2_h12 = conv(x_s_noise_fading2,g12);

%% calculate the channel information of the virtual end-to-end link
% firstly, we obtain the 10 frames CP data
x_s_noise_fading1_cp = x_s_noise_fading1(9:32);
for n = 1:1:9
    x_s_noise_fading1_cp =[x_s_noise_fading1_cp(:)' x_s_noise_fading1(160*n+9:160*n+32)']';
end

x_s_noise_fading1_h21_cp = x_s_noise_fading1_h21(9:32);
for n = 1:1:9
    x_s_noise_fading1_h21_cp =[x_s_noise_fading1_h21_cp(:)' x_s_noise_fading1_h21(160*n+9:160*n+32)']';
end

x_s_noise_fading2_cp = x_s_noise_fading2(9:32);
for n = 1:1:9
    x_s_noise_fading2_cp =[x_s_noise_fading2_cp(:)' x_s_noise_fading2(160*n+9:160*n+32)']';
end

x_s_noise_fading2_h12_cp = x_s_noise_fading2_h12(9:32);
for n = 1:1:9
    x_s_noise_fading2_h12_cp =[x_s_noise_fading2_h12_cp(:)' x_s_noise_fading2_h12(160*n+9:160*n+32)']';
end

v12 = x_s_noise_fading1_cp.*x_s_noise_fading2_h12_cp;
v21 = x_s_noise_fading2_cp.*x_s_noise_fading1_h21_cp;
% important record: the square of the mean is not equal to the mean of the
% square as follown
%v12_mul_pwr = mean(abs(x_s_noise_fading1_cp).^2)*mean(abs(x_s_noise_fading2_h12_cp).^2);

v12_mul_pwr = mean((abs(x_s_noise_fading1_cp).^2).*(abs(x_s_noise_fading2_h12_cp).^2));
v21_mul_pwr = mean((abs(x_s_noise_fading2_cp).^2).*(abs(x_s_noise_fading1_h21_cp).^2));

v12_pwr = mean(abs(v12).^2);
v21_pwr = mean(abs(v21).^2);