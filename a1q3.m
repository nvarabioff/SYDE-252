% SYDE 252 %
% MATLAB Assignment 1 %

% Problem 3: Numeric Approximation of Fourier Transform of Continuous Time Signals

% * Make helper function for finding period for t range and finding
% possible freq's with t value, and vice versa for ift
% * Better way to validate test cases
% * More test cases
% * Compare test cases with built in Fourier Transform functions

clc, clear;
clf;

Fs = 1000; %Sampling frequency
Ts = 1/Fs; % Sampling Period
L = 1000; % Length of signal
t = (0:L-1)*Ts; % Time vector
L2 = 2^nextpow2(L);
w = (0:(Fs/L2):(Fs/2-Fs/L2));
dim = 2;

% Test 1: sin(100t)
x = sin(100*t);

% Test 2: complex exponential
% x = exp(1i*80*t);

% Test 3: cos(60t)
% x = cos(60*t);

xw = MyFT(x,t,w);
xt = MyiFT(xw,w,t);

figure(1);
subplot(5,1,1);
hold on;
title('x(t)');
xlabel('t');
plot(t, x);
hold off;

subplot(5,1,2);
hold on;
mag1 = (sqrt(real(xw).^2 + imag(xw).^2));
% plot(w,real(xw));
% plot(w,imag(xw));
plot(w, mag1);
% title('Fourier Transform of x(t) to X(w): Blue Real, Red Imaginary');
title('Fourier Transform of x(t) to X(w)');
xlabel('w(rad/s)');
hold off;

subplot(5,1,3);
hold on;
Y = fft(x, L2, dim);
P2 = abs(Y);
P1 = P2(:,1:L2/2+1);
plot(w, P1(1: L2/2));
title('Built-in Fourier Transform of x(t) to X(w)');
xlabel('f(Hz)');
hold off;

subplot(5,1,4);
hold on;
%mag2 = sqrt(real(xt).^2 + imag(xt).^2);
plot(t, real(xt)/1000);
%plot(t,imag(xt));
%plot(t, mag2);
title('Inverse Fourier Transform of X(w) back to x(t): Blue Real, Red Imaginary');
xlabel('t');
hold off;

% subplot(5,1,5);
% hold on;
% X = ifft(Y, L2, dim);
% plot(t, X);
% title('Built-in Inverse Fourier Transform of X(w) back to x(t)');
% xlabel('t');
% hold off;

% Fourier Transform (Analysis Equation)
% Parameters:
    % Xt is the function in time domain
    % t is sampling rate and range of CT signal in time domain.
        % Rate should be high (i.e. small intervals = 0.01)
        % Range must be at least range of period, but should be higher
    % w is range of frequencies and display rate
        % Rate should be same as t (i.e. small intervals = 0.01)
        % Range should be t/2*pi
function [Xw] = MyFT(Xt, t, w)
    
    Xw = zeros(1, length(w));   % make Xw with range encompassing possible frequencies
    
    for i = 1:length(w) % for every value on frequency representation Xw
        v = 0;  % v is value at current index of Xw
        for j = 1:length(t) % sum every value of time representation * comlpex exponential
            v = v + Xt(j)*exp(-1i*t(j)*w(i));
        end
        Xw(i) = v;  % assign calculated value to index of Xw
    end
end

% Inverse Fourier Transform (Synthesis Equation)
    % Parameters:
    % Xw is the function in frequency domain.
    % w is range of frequencies sampling rate of fourier transform Xw
        % Rate should be high (i.e. small intervals = 0.01)
        % Range should be all possible frequencies of signal
    % t is range of signal displayed in time domain and 
        % Rate should be same as w (i.e. small intervals = 0.01)
        % Range must be 2*pi*w
% Inverse Fourier Transform is same operations as Fourier Transform other
% than the conjugate of the complex exponential, the division by 2pi, and
% the integrating variable. Therefore, we can use the same function
% structure with a few minor changes.
function [Xt] =  MyiFT(Xw, w, t)
  
  Xt = zeros(1,length(t));  % representation in time domain must be same length as possible periods
  
  for i = 1:length(t)   % Switched looping variable because different integrating variable in inverse fourier transforms
    v = 0;
    for j = 1:length(w)
      v = v + Xw(j)*exp(1i*w(j)*t(i));
    end
    Xt(i) = v/(2*pi);   % Remember division by 2pi
  end
end
