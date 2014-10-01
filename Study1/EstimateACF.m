function [ ryp ] = EstimateACF(x,type)
%Takes the signal x and returns the estimated ACF ryp at k
%Estimation by Blackman-Tukey's  and Bartlett's method 
N = max(size(x));
ryp = zeros(1,N);

for k = -N/2+1:N/2
    s = 0;
    for n = 1:(N-abs(k))
        s = s+x(n+abs(k))*x(n);
    end
    ryp(k+N/2) = s;
    switch type
    case 'BmanT'
        ryp(k+N/2) = ryp(k+N/2)/(N-abs(k));
    case 'Blett'
        ryp(k+N/2) = ryp(k+N/2)/(N); 
    end
end     
end

