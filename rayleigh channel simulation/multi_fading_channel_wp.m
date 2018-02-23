

sampleRate500KHz = 500e3;    	% Sample rate of 500K Hz
sampleRate20KHz  = 20e3;     	% Sample rate of 20K Hz
maxDopplerShift  = 200;      	% Maximum Doppler shift of diffuse components (Hz)
delayVector = (0:5:15)*1e-6; 	% Discrete delays of four-path channel (s)
gainVector  = [0 -3 -6 -9];  	% Average path gains (dB)

KFactor = 10;            		% Linear ratio of specular power to diffuse power
specDopplerShift = 100;  		% Doppler shift of specular component (Hz)

% Configure a Rayleigh channel object
rayChan = comm.RayleighChannel( ...
    'SampleRate',          sampleRate500KHz, ...
    'PathDelays',          delayVector, ...
    'AveragePathGains',    gainVector, ...
    'MaximumDopplerShift', maxDopplerShift, ...
    'RandomStream',        'mt19937ar with seed', ...
    'Seed',                10, ...
    'PathGainsOutputPort', true);

% Configure a Rician channel object
ricChan = comm.RicianChannel( ...
    'SampleRate',              sampleRate500KHz, ...
    'PathDelays',              delayVector, ...
    'AveragePathGains',        gainVector, ...
    'KFactor',                 KFactor, ...
    'DirectPathDopplerShift',  specDopplerShift, ...
    'MaximumDopplerShift',     maxDopplerShift, ...
    'RandomStream',            'mt19937ar with seed', ...
    'Seed',                    100, ...
    'PathGainsOutputPort',     true);

qpskMod = comm.QPSKModulator( ...
    'BitInput',    true, ...
    'PhaseOffset', pi/4);

% Number of bits transmitted per frame is set to be 1000. For QPSK
% modulation, this corresponds to 500 symbols per frame.
bitsPerFrame = 1000;
msg = randi([0 1],bitsPerFrame,1);

% Modulate data for transmission over channel
modSignal = qpskMod(msg);

% Apply Rayleigh or Rician channel object on the modulated data
rayChan(modSignal);
ricChan(modSignal);

release(rayChan);
release(ricChan);

rayChan.Visualization = 'Impulse and frequency responses';
rayChan.SamplesToDisplay = '100%';

numFrames = 2;

for i = 1:numFrames % Display impulse and frequency responses for 2 frames
    % Create random data
    msg = randi([0 1],bitsPerFrame,1);
    % Modulate data
    modSignal = qpskMod(msg);
    % Filter data through channel and show channel responses
    rayChan(modSignal);
end


release(rayChan);
rayChan.Visualization = 'Doppler spectrum';

numFrames = 5000;

for i = 1:numFrames % Display Doppler spectrum from 5000 frame transmission
    msg = randi([0 1],bitsPerFrame,1);
    modSignal = qpskMod(msg);
    rayChan(modSignal);
end

release(rayChan);
rayChan.Visualization = 'Impulse and frequency responses';
rayChan.SampleRate = sampleRate20KHz;
rayChan.SamplesToDisplay = '25%';  % Display one of every four samples

numFrames = 2;

for i = 1:numFrames % Display impulse and frequency responses for 2 frames
    msg = randi([0 1],bitsPerFrame,1);
    modSignal = qpskMod(msg);
    rayChan(modSignal);
end


release(rayChan);
rayChan.PathDelays = 0;        % Single fading path with zero delay
rayChan.AveragePathGains = 0;  % Average path gain of 1 (0 dB)

for i = 1:numFrames % Display impulse and frequency responses for 2 frames
    msg = randi([0 1],bitsPerFrame,1);
    modSignal = qpskMod(msg);
    rayChan(modSignal);
end

release(rayChan);
rayChan.Visualization = 'Off'; % Turn off System object's visualization
ricChan.Visualization = 'Off'; % Turn off System object's visualization

% Same sample rate and delay profile for the Rayleigh and Rician objects
ricChan.SampleRate       = rayChan.SampleRate;
ricChan.PathDelays       = rayChan.PathDelays;
ricChan.AveragePathGains = rayChan.AveragePathGains;

% Configure a Time Scope System object to show path gain magnitude
gainScope = dsp.TimeScope( ...
    'SampleRate', rayChan.SampleRate, ...
    'TimeSpan',   bitsPerFrame/2/rayChan.SampleRate, ... % One frame span
    'Name',       'Multipath Gain', ...
    'ShowGrid',   true, ...
    'YLimits',    [-40 10], ...
    'YLabel',     'Gain (dB)');

% Compare the path gain outputs from both objects for one frame
msg = randi([0 1],bitsPerFrame,1);
modSignal = qpskMod(msg);
[~, rayPathGain] = rayChan(modSignal);
[~, ricPathGain] = ricChan(modSignal);
% Form the path gains as a two-channel input to the time scope
gainScope(10*log10(abs([rayPathGain, ricPathGain]).^2));

clear hRicChan hMultipathGain;
release(rayChan);

rayChan.PathDelays = delayVector;
rayChan.AveragePathGains = gainVector;
rayChan.MaximumDopplerShift = 5;

% Configure a Constellation Diagram System object to show received signal
constDiag = comm.ConstellationDiagram( ...
    'Name', 'Received Signal After Rayleigh Fading');

numFrames = 16;

for n = 1:numFrames
    msg = randi([0 1],bitsPerFrame,1);
    modSignal = qpskMod(msg);
    rayChanOut = rayChan(modSignal);
    % Display constellation diagram for Rayleigh channel output
    constDiag(rayChanOut);
end


release(rayChan);
release(constDiag);
rayChan.SampleRate = sampleRate500KHz;

for n = 1:numFrames
    msg = randi([0 1],bitsPerFrame,1);
    modSignal = qpskMod(msg);
    rayChanOut = rayChan(modSignal);
    constDiag(rayChanOut);
end