%% this example process a binary data stream using a communication system that consists of 
%  a baseband modulator, channel, and demodulator. 
%% The system's bit error rate (BER) is computed and the transmitted and received signals 
%  are displayed in a constellation diagram.
% ********************************************

% the qammod function dose not apply pulse shaping
clc;
clear;
close all;
M = 16; 		% Size of signal constellation
k = log2(M);	% number of bits per symbol
n = 3e5;		% numberi of bits to process
numSamplePerSymbol = 1; % oversampling factor

rng default;    % use default random number generator
dataIn = randi([0 1],n,1); %generate vector of binary data

% plot the first 40 bits in a stem plot
stem(1:40,dataIn(1:40),'filled');
title('random bits');
xlabel('Bit index');
ylabel('Binary value');


% convert teh binary signal to an integer-valued signal
dataInMatrix = reshape(dataIn,length(dataIn)/k,k);
dataSymbolsIn = bi2de(dataInMatrix);

% plot the first 10 symbols
stem(dataSymbolsIn(1:10));
title('Random Symbols');
xlabel('Symbol Index');
ylabel('Integer Value');

% Modulate using 16-QAM
dataMod = qammod(dataSymbolsIn,M,'bin');		%binary coding, phase offset = 0
dataModG = qammod(dataSymbolsIn,M);				% Gray coding, Phase offset = 0;

%sPlotFig = scatterplot(dataMod,1,0,'r.');

% calculate the SNR when the channel has an Eb/N0 = 10 dB.
EbNo = 10;
snr = EbNo + 10*log10(k) - 10*log10(numSamplePerSymbol);

% pass the signal through AWGN channel
receivedSignal = awgn(dataMod,snr,'measured');
receivedSignalG = awgn(dataModG,snr,'measured');

% create the constellation diagram
sPlotFig = scatterplot(receivedSignal,1,0,'g.');
hold on;
scatterplot(dataMod,1,0,'k*',sPlotFig);

% demodulate 16-QAM
dataSymbolsOut = qamdemod(receivedSignal,M,'bin');
dataSymbolsOutG = qamdemod(receivedSignalG,M);

% convert the integer-valued signal to binary signal
dataOutMatrix = de2bi(dataSymbolsOut,k);
dataOut = dataOutMatrix(:);
dataOutMatrixG = de2bi(dataSymbolsOutG,k);
dataOutG = dataOutMatrixG(:);

% compute teh system BER
[numErrors,ber] = biterr(dataIn,dataOut);
fprintf('\nThe binary coding bit error rate = %5.2e, base on %d errors(total bits %d)\n', ber, numErrors,n);

[numErrorsG,berG] = biterr(dataIn,dataOutG);
fprintf('\nThe binary coding bit error rate = %5.2e, base on %d errors(total bits %d)\n', berG, numErrorsG,n);


% plot signal constellations 画出信号调制星座图
M = 16;                         % Modulation order
x = (0:15)';                    % Integer input

y1 = qammod(x,16,'bin');
scatterplot(y1)
text(real(y1)+0.1, imag(y1), dec2bin(x))
title('16-QAM, Binary Symbol Mapping')
axis([-4 4 -4 4])

y2 = qammod(x,16,'gray');  % 16-QAM output, Gray-coded
scatterplot(y2)
text(real(y2)+0.1, imag(y2), dec2bin(x))
title('16-QAM, Gray-coded Symbol Mapping')
axis([-4 4 -4 4])

% pulse shaping using a rasied consine filter (a pair of square-root raised consine (RRC) filter)
% 平方根升余弦滤波器（RRC）
% 用来做signal shaping的，目的在一定的带宽要求下，尽量减少码间串扰（ISI）
% 匹配滤波的目标也是为了修正ISI带来的信号畸变。
% 升余弦滚降信号用来消除码间串扰，实际实现时采用的方式是由发送端的基带成形滤波器和接收端的匹配滤波器
% 两个环节共同实现。传输系统的传递函数为二者的乘积，所以每个环节均为平方根升余弦滤波器，这样可以降低滤波器的是实现难度

% 在数字通信中，实际发射出的信号是各个离散样值序列通过成形滤波器后的成形脉冲序列，匹配滤波器是为了使得
% 抽样时刻信噪比最大。
% 在发送端成形滤波器是根余弦滤波器，接收端同样使用根余弦滤波器匹配滤波时，既能够使得抽样时刻信噪比
% 最高（即完成匹配滤波器的作用），又能够在一定的带限平坦信道中不引入码间干扰（满足Nyquist无码间干扰准则）


% establis simulation framework
M = 16;							% size of signal constellation
k = log2(M);					% number of bits per symbol
numBits = 3e5;
numSamplePerSymbol = 4

% create raised consine filter
span = 10;						% filter span in symbols
rolloff = 0.25;					% rolloff factor of filter

rrcFilter = rcosdesign(rolloff,span,numSamplePerSymbol);
fvtool(rrcFilter,'Analysis','Impulse');	% display teh RRC filter response

% BER simulation
rng default;					% use default random number generator
dataIn = randi([0 1], numBits, 1); % generate vector of binary data

dataInMatrix = reshape(dataIn,length(dataIn)/k,k);
dataSymbolsIn = bi2de(dataInMatrix);

dataMod = qammod(dataSymbolsIn,M);

% 滤波器的过采样和欠采样,
txSignal = upfirdn(dataMod,rrcFilter,numSamplePerSymbol,1);

EbNo = 10;
snr = EbNo + 10*log10(k) - 10*log10(numSamplePerSymbol);

% pass the signal through an AWGN channel
rxSignal = awgn(txSignal,snr,'measured');

rxFiltSignal = upfirdn(rxSignal,rrcFilter,1,numSamplePerSymbol);
rxFiltSignal = rxFiltSignal(span+1:end-span);

dataSymbolsOut = qamdemod(rxFiltSignal,M);

dataOutMatrix = de2bi(dataSymbolsOut,k);
dataOut = dataOutMatrix(:);

% 从结果分析，错误比特差不多
[numErrors, ber] = biterr(dataIn, dataOut);
fprintf('\nThe bit error rate = %5.2e, based on %d errors\n', ...
    ber, numErrors);

% create a scatter plot of the received signal before and after filtering
h = scatterplot(sqrt(numSamplePerSymbol)*...
    rxSignal(1:numSamplePerSymbol*5e3),...
    numSamplePerSymbol,0,'g.');
hold on;
scatterplot(rxFiltSignal(1:5e3),1,0,'kx',h);
title('Received Signal, Before and After Filtering');
legend('Before Filtering','After Filtering');
axis([-5 5 -5 5]); % Set axis ranges
hold off;



%% error correction using a convolutional code (纠错)
clc;
clear;
close all;
% establish simulation framework
M = 16;
k = log2(M);
numBits = 1e5;				% 比特数，每四个比特为一个码
numSamplePerSymbol = 4; 	% 每个码采样率

rng default;
dataIn = randi([0 1], numBits,1);

% define a convolutional coding trellis for a rate 2/3 code
tPoly = poly2trellis([5 4],[23 35 0; 0 5 13]);
codeRate = 2/3;
dataEnc = convenc(dataIn,tPoly);

dataEncMatrix = reshape(dataEnc,length(dataEnc)/k,k);
dataSymbolsIn = bi2de(dataEncMatrix);

dataMod = qammod(dataSymbolsIn,M);

% create raised consine filter
span = 10;						% filter span in symbols
rolloff = 0.25;					% rolloff factor of filter

rrcFilter = rcosdesign(rolloff,span,numSamplePerSymbol);

txSignal = upfirdn(dataMod,rrcFilter,numSamplePerSymbol,1);

EbNo = 10;
snr = EbNo + 10*log10(k) - 10*log10(numSamplePerSymbol);

% pass the signal through an AWGN channel
rxSignal = awgn(txSignal,snr,'measured');

rxFiltSignal = upfirdn(rxSignal,rrcFilter,1,numSamplePerSymbol);
rxFiltSignal = rxFiltSignal(span+1:end-span);

dataSymbolsOut = qamdemod(rxFiltSignal,M);

dataOutMatrix = de2bi(dataSymbolsOut,k);
codedDataOut = dataOutMatrix(:);
traceBack = 16;
numCodeWords = floor(length(codedDataOut)*2/3);
dataOut = vitdec(codedDataOut(1:numCodeWords*3/2),tPoly,traceBack,'cont','hard');

decDelay = 2*traceBack;                                     % Decoder delay, in bits
[numErrors, ber] = ...
   biterr(dataIn(1:end-decDelay),dataOut(decDelay+1:end));       

fprintf('\nThe bit error rate = %5.2e, based on %d errors\n', ...
    ber, numErrors)