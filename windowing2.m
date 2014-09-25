function [y] = windowing2(acf, P)
%WINDOWING2 Summary of this function goes here
%   Detailed explanation goes here
N = length(acf);

w2 = window(@blackmanharris,P);
%w2 = window(@triang,P);

w = zeros(1,N);
w(ceil((N-P)/2):floor((N+P)/2)) = w2;

y = abs(fft(acf.*w));


end

