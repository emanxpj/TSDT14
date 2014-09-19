function [ ryp ] = EstimateACF(x,type)
%Takes the signal x and returns the estimated ACF ryp at k
%Estimation by Blackman-Tukey's  and Bartlett's method 
N = max(size(x));
ryp = zeros(1,2*N);

for k = -N+1:N
    s = 0;
    for n = 1:(N-abs(k))
        s = s+x(n+abs(k))*x(n);
    end
    ryp(k+N) = s;
    switch type
    case 'BmanT'
        ryp(k+N) = ryp(k+N)/(N-abs(k));
    case 'Blett'
        ryp(k+N) = ryp(k+N)/(N); 
    end
end     
end

