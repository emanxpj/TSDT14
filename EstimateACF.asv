function [ ryp ] = EstimateACF(x,k,type)
%Takes the signal x and returns the estimated ACF ryp at k
%Estimation by Blackman-Tukey's  and Bartlett's method 
N = max(size(x));
ryp = zeros(1,N);
s = 0;
for k = 1:N
    for n = 1:(N-abs(k))
        s = s+x(n+abs(k))*x(n);
    end
    ryp(k) = s;
    switch type
    case 'BmanT'
        ryp(k) = ryp(k)/(N-abs(k));
    case 'Blett'
        ryp(k) = ryp(k)/N; 
    end
end     
end

