% SYDE 252 %
% MATLAB Assignment 1 %

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