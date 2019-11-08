% Problem 1: Numeric approximation of continuous time convolution


% Problem 3: Numeric Approximation of Fourier transform of continuous time signals

%% Fourier Transform (Analysis Equation)
  % Parameters:
    % Xt = function
    % t = 
    % w = 
  % Returns:
    % Xw = the Fourier transform
function [Xw] = MyFT(Xt, t, N)
  % Define complex exponential. i is complex number, a is t, b is N
  e = @(a,b) exp(-1i*a*b);
  
  Xw = zeros(1, length(N));
  
  % First loop represents the summation from 0 to N - 1
  for i = 1:length(N)
    v = 0; 
    % Second loop represents the value at each non-zero x(t) * e^-jwn
    for j = 1:length(t)
      v = v + Xt(j)*e(t(j),N(i));
    end
    Xw(i) = v;
  end
end

%% Inverse Fourier Transform (Synthesis Equation)

function [Xt] =  MyiFT(Xw, w, t)
  % Define complex exponential. i is complex number, a is time period t, b is frequency w
  e = @(a,b) exp(1i*a*b);
  
  Xt = zeros(1,length(t));
  
  for i = 1:length(t)
    v = 0;
    for j = 1:length(w)
      v = v + Xw(j)*e(w(j),t(i));
    end
    Xt(i) = v/(2*pi);
  end
end





