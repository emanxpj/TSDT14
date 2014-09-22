function [ PerAv ] = PerAv(x,k)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
y = Periodogram(x);

N = max(size(x));
l = N/k;
PerAv = zeros(1,N);

for p = 1:k
    PerAv((p-1)*l+1:p*l) = mean(y(l*(p-1)+1:p*l));



end

