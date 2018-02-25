% AM modulation

fs = 1e5;
t = (0:1/fs:100)';

fc = 1e3;
x = sin(2*pi*t);

ydouble = ammod(x,fc,fs);
ysingle = ssbmod(x,fc,fs);

sa = dsp.SpectrumAnalyzer('SampleRate',fs, ...
    'PlotAsTwoSidedSpectrum',false, ...
    'YLimits',[-60 40]);
step(sa,ydouble)

step(sa,ysingle)