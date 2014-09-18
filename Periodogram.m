function [ Ryp ] = Periodogram( x,k )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
N = max(size(x));
s = 0;
Ryp = zeros(1,N);
X = zeros(1,N);
for k = 1:N
    for n = 1:N
        s =s + x(n)*exp(-i*2*pi*k*n);
    end
    X(k) = s;
end
Ryp = (1/N)*abs(X).^2;
end

