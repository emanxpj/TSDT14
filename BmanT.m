function [ ryp ] = BmanT(x,k,dt)
%Takes the signal x and returns the estimated ACF ryp at k with sample rate
%dt.
%   Estimation by Blackman-Tukey's method 

N = size(x);
s = 0;
if(abs(k) < N(1))
    for n = 1:(N(1)-abs(k)/dt-1)
        s = s+x(logical(n+abs(k/dt)))*x(n);
    end
    ryp = s;
else
    ryp = 0;
end
end

