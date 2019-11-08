% Problem 1: Numeric Approximation of continuous-time convolution %

clc, clear;
clf;

% Testing

% Test 1:
% x = ones(1,100);
% y = 0.01:0.01:2;

% % Test 2:
% t = 0:0.01:8*pi;
% x = sin(t);
% y = ones(1,10);

% Test 3:
t = 0:0.01:8*pi;
x = sin(t);
y = -2 * ones(1,10);

% Testing against built-in MATLAB function
convSoln = conv(x,y);
myConvSoln = myConvolution(x,y);
% disp("Using Conv function");
% disp(convSoln);

% Plot
figure(1);
subplot(4,1,1);
hold on;
title('x(t)');
xlabel('t');
plot(x);
hold off;

subplot(4,1,2);
hold on;
title('y(t)');
xlabel('t');
plot(y);
hold off;

subplot(4,1,3);
hold on;
title('My Convolution');
xlabel('t');
plot(myConvSoln);
hold off;

subplot(4,1,4);
hold on;
title('Built in convolution');
xlabel('t');
plot(convSoln);
hold off;

function [z] = myConvolution(x,y)
    
    % Create an array of zeros "z" to represent the solution z(t)
    % The array z is the length of x plus the length of y - 1
    % This is because the value of z(t) is only nonzero for a period of t
    % while x and y overlap. This period will only have a range of (length 
    % of y) + (length of x), because shifting further than this in either 
    % direction will result in z having a value of zero, because either x 
    % or y will be zero. Therefore, we only need to define the solution 
    % array z from 1 to (length of x) + (length of y) - 1.
    z = zeros(1, length(x)+length(y)-1);
    
    % Iterate through the array z to store the solution at that index
    % This index defines how far shifted the array y is
    for t = 1:length(z)
        v = 0;  %this will be the value at the index of the result array z
        for k = 1:t
            % Skip this iteration if t is out of range of length of x (i.e. x at that index = 0)
            if k > length(x)
                continue;
            % Skip this iteration if t - k + 1 is out of range of length of y (i.e. y at that index = 0)
            elseif t - k + 1 > length(y)    
                continue;
            end
            % Set the value of v (the value at that index of the solution to be x(t)*y(t - k + 1)
            v = v + x(k)*y(t - k + 1);  
        end
        % Assign the value calculated by the above loops to the index of the solution array
        z(t) = v;
    end
    
end


% Problem 3: Numeric Approximation of Fourier Transform of Continuous Time
% Signals

% * Make helper function for finding period for t range and finding
% possible freq's with t value, and vice versa for ift
% * Better way to validate test cases
% * More test cases

clc, clear;
clf;

% Test 1: sin(t)
% t = 0:0.01:8*pi;
% w = 0:0.01:4;
% x = sin(t);

% Test 2: cos(3t)
% t = 0:0.1:24*pi;
% w = 0:0.01:12;
% x = cos(3*t);

% Test 3: complex exponential
t = 0:0.01:8*pi;
w = 0:0.01:4;
x = exp(i*2*t);

xw = MyFT(x,t,w);
xt = MyiFT(xw,w,t);

figure(2);
subplot(3,1,1);
hold on;
title('x(t)');
xlabel('t');
plot(t,x);
hold off;

subplot(3,1,2);
hold on;
mag1 = sqrt(real(xw).^2 + imag(xw).^2);
% plot(w,real(xw));
% plot(w,imag(xw));
plot(w, mag1);
% title('Fourier Transform of x(t) to X(w): Blue Real, Red Imaginary');
title('Fourier Transform of x(t) to X(w): Magnitude');
xlabel('w');
hold off;

subplot(3,1,3);
hold on;
mag2 = sqrt(real(xt).^2 + imag(xt).^2);
plot(t,real(xt));
plot(t,imag(xt));
% plot(t, mag2);
title('Inverse Fourier Transform of X(w) back to x(t): Blue Real, Red Imaginary');
xlabel('t');
hold off;

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




