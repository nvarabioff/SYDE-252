% SYDE 252 %
% MATLAB Assignment 1 %

% Problem 4 %
% Option 2: Identification of Piano Keys %

% Differences bewteen Built in and Custom:
    % fft_length is multiplied by 2pi on line 91 and 99
    % Fourier Transform line 91 - no need to split Xjw because custom function does not mirror impulse response

% Toggle music 1 and music 2 with line 20

clc, clear;
clf;

% Import sound file
load('pianoNotes.mat');
load('Music1.mat'); % Fur Elise
music1Data = acqData;   % Variable storing music1 data
load('Music2.mat'); % Jingle Bells
music2Data = acqData;   % Variable storing music 2 data
currData = music1Data;  % Toggle this to test Music1 and Music2 files

% Convert note series to function Xt
xt_data = currData(:,1);
xt_length = length(xt_data);
Xt = zeros(1,xt_length);
for i = 1:xt_length
    Xt(i) = currData(i);
end

% Plot function Xt
figure(1);
hold on;
title('x(t)');
xlabel('t');
plot(Xt);
hold off;

% Information for Fourier Transform
Fs = 16000; % Sampling frequency
Ts = 1/Fs;  % Sampling Period
xt_length2 = 2^nextpow2(xt_length); % Creates new input length that is next power 2 of original length
fft_length = Fs*(0:xt_length2/2-1)/xt_length2;  % Length of signal in frequency domain

% Find note lengths
noteLengths = [];   % Array storing the lengths of each note
lastNote = 0;   % Index of last note
nLength = Fs/4; % Length of note. Stop detecting notes for 1/4 second after note detected to remove noise
for k = 1:xt_length-1
    if lastNote ~= 0    % Stop detecting notes for 1/4 second after note detected to remove noise
        if k - lastNote < nLength
            continue;
        end
    end
    if abs(Xt(k+1) - Xt(k)) > 0.03     % Find peaks (i.e. Notes)
        noteLengths = [noteLengths; k - lastNote];  % Add the length of the note to array
        lastNote = k;   % Store index of current note
    end
end

% noteLengths = [noteLengths(2:end); xt_length-lastNote-5000];  % Remove first "note" (noise) and add last note. Need the -5000 or the last note will have the wrong frequency
noteLengths = [noteLengths(2:end); 8000];

% Fourier Transform and detect frequency of each note
currentNote = [];   % Current note
noteLength = 0; % Length of each note
counter = 1;    % Number of current note
noteFrequencies = [];   % Array of note frequencies
recording = false;
for j = 1:xt_length-1   % Loop through entire data
    if recording
        currentNote = [currentNote; Xt(j)];
        noteLength = noteLength - 1;
        if noteLength > 0
            continue;
        elseif noteLength <= 0
            recording = false;
            
            t = 0:Ts:(length(currentNote)-1)/Fs;  % Length of time domain vector for the current note
            % Plot the current note
            figure(counter);
            subplot(2,1,1);
            hold on;
            title('x(t)');
            xlabel('t');
            plot(t,currentNote);
            hold off;
            
            % Custom Fourier Transform
            Xjw = MyFT(currentNote, t, fft_length*2*pi);
            
            % Plot Fourier transform
            figure(counter);
            subplot(2,1,2);
            hold on;
            title('X(jw)');
            xlabel('f');
            plot(fft_length*2*pi,abs(Xjw));  % Custom
            hold off;
            
            % Find the maximum of the Fourier Transform - first value above x
            m = 1;
            while abs(Xjw(m+1) - Xjw(m)) < 2.63   % Detects first big spike in frequency (Built In)
                m = m + 1;  % Index of the spike is the frequency
            end
            
            noteFrequencies = [noteFrequencies; m/8.192];    % Add to noteFrequencies list, normalized
            
            currentNote = [];
        end
    elseif abs(Xt(j+1) - Xt(j)) > 0.03
        recording = true;
        noteLength = noteLengths(counter)-1;
        counter = counter + 1;
        continue;
    end
end
            
% For each note frequency, find the note and add it to the array of note letters
noteLetters = [];
for p = 1:length(noteFrequencies)
    noteLetters = [noteLetters; findNote(noteFrequencies(p), noteFreqs, noteNamesFull)];
end

% My Fourier Transform Function
function [Xw] = MyFT(Xt, t, w)
    Xw = zeros(1, length(w));   % make Xw with range encompassing possible frequencies
    for i = 1:length(w) % for every value on frequency representation Xw
        v = 0;  % v is value at current index of Xw
        for j = 1:length(t) % sum every value of time representation * comlpex exponential
            v = v + Xt(j)*exp(-1i*t(j)*w(i));
        end
        Xw(i) = v;  % assign calculated value to index of Xw. Added the /2pi cause of radians to Hz
    end
end

% Function to find notes from given data
function [noteLetter] = findNote(freq, allNotes, noteNames) % allNotes is the given note frequencies, noteNames is the note names
    for i = 2:length(allNotes)-1
        if freq < allNotes(i)   % Find closest note from given data
            if abs(freq - allNotes(i-1)) < abs(freq - allNotes(i))
                noteLetter = noteNames(i-1);
            else
                noteLetter = noteNames(i);
            end
            break;
        end
    end
end

    
