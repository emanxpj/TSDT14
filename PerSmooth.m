function [ y2 ] = PerSmooth(x,N)
%PERSMOOTH Summary of this function goes here
%   Detailed explanation goes here

y = Periodogram(x);
k = max(size(x))
w = linspace(-k/2,k/2,k);
S = 1/N*sinc(N*w).^2;
% assign the wn vector as a triangle then transform?
y2 = cconv(S,y,2*k);


end

