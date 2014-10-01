function [y2] = PerSmooth2(x,N)
%PERSMOOTH2 Summary of this function goes here
%   Detailed explanation goes here

k = max(size(x)); 
w = linspace(-1/2,1/2,k);
S = sinc(w*N).^2;
plot(w,S);

y = movingAverage(Periodogram(x),5);
y2 = movingAverage(y,5);

end

function y = movingAverage(x, N)
   n = linspace(-5,5,11);
   k = zeros(1,11);
   k(abs(n) < 5) = (1-abs(n(abs(n)<5)/5))/5;
   
   y = cconv(x,k,max(size(x)));
   %k = ones(1, N) / N;
   %y = cconv(x, k, max(size(x)));
end

function y = movingAverage3(x, N)
   k = ones(1, N) / N;
   y = cconv(x, k, max(size(x)));
end

function y = movingAverage2(x, N)
    %No idea why this doesn't work. 
    
    w = linspace(-1/2,1/2,max(size(x)));
    S = 1/N*sinc(w*N);
    y = cconv(x, S, max(size(x)));
end
