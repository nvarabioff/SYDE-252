load('Music1.mat');
filename = 'Music1.wav';
y = acqData;
Fs = 22100; % This value is important
audiowrite(filename, y, Fs);
clear y Fs;
