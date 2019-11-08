% Problem 1: Numeric approximation of continuous time convolution


% Problem 3: Numeric Approximation of Fourier transform of continuous time signals

%% Fourier Transform
  % Parameters:
  % Xt = function
  % t = fundamental period -> n in discrete
  % w = fundamental frequency
function [Xw] = MyFT(Xt, t, w)
  % Define complex exponential. i is iterator, a is time period t, b is frequency w
  e = @(a,b) exp(-1i*a*b);
  
  Xw = zeros(1, length(w));
  for i = 1:length(w)
    v = 0;
    for j = 1:length(t)
      v = v + Xt(j)*e(t(j), w(i));
    end
    Xw(i) = v;
  end
end

%% Inverse Fourier Transform
