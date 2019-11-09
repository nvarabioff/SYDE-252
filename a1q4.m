% SYDE 252 %
% MATLAB Assignment 1 %

% Problem 4 %
% Option 2: Identification of Piano Keys %

clc, clear;
clf;

% Import sound file

load('pianoNotes.mat');
load('Music1.mat');
music1Data = acqData;
load('Music2.mat');
music2Data = acqData;

currData = music2Data;

% Convert note series to function Xt

xt_data = currData(:,1);
xt_length = length(xt_data);
Xt = zeros(1,xt_length);

for i = 1:xt_length
    Xt(i) = currData(i);
end

% Plot function Xt
figure(1);
subplot(2,1,1);
hold on;
title('x(t)');
xlabel('t');
plot(Xt);
hold off;

% Detect note, split notes

currentNote = [];   % Current note
noteLength = 0; % Standard length of each note
counter = 1;    % Number of current note
for j = 1:xt_length-1   % Loop through entire data
    if noteLength > 0  % If recording a note
        currentNote = [currentNote; Xt(j)]; % Dynamically add to the current note for 5000 units
        noteLength = noteLength - 1;
    end
    if Xt(j) > 0.025 && noteLength == 0 % Detecting the start of a new note
        noteLength = 5000;
        counter = counter + 1;
        
        if counter > 2
            % Plot the current note
            figure(counter);
            subplot(2,1,1);
            hold on;
            title('x(t)');
            xlabel('t');
            plot(currentNote);
            hold off;

            % Fourier transform
            Xjw = fft(currentNote);

            % Plot Fourier transform
            figure(counter);
            subplot(2,1,2);
            hold on;
            title('X(jw)');
            xlabel('f');
            plot(abs(Xjw));
            hold off;
        end
        
        currentNote = [];   % Clear the value of current note
    end
end

% Plot Notes
figure(1);
subplot(2,1,2);
hold on;
title('x(t)');
xlabel('t');
plot(Xt(2));
hold off;

% Fourier Transform each note
% Graph the Fourier Transform
% Detect the frequency of the current note
% Match it with a note in the pianoNotes.mat file
% Display the note name
