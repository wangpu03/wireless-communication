%% this example process a binary data stream using a communication system that consists of 
%  a baseband modulator, channel, and demodulator. 
%% The system's bit error rate (BER) is computed and the transmitted and received signals 
%  are displayed in a constellation diagram.
% ********************************************

clc;
clear;
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
