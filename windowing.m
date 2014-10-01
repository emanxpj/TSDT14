function [y] = windowing(x,n,method)
%WINDOWING Summary of this function goes here
%   Detailed explanation goes here

switch method
    case 'square'
        w = ones(1,n)/n;
        x1 = Periodogram(x);
        y = filter(w,1,x1);
    case 'hamming'
        w = hanning(n)/n;
        x1 = Periodogram(x);
        y = filter(w,1,x1);    
    case 'pure'
        y = PeriodFourier(x)
        
    case 'mat'
        w = linspace(-1/2,1/2,length(x));
        y  = smooth(w,Periodogram(x),0.1);
        
    case 'mat2'
        width = n; 
        x = Periodogram(x);
        y = filter(ones(width,1)/width,1,x);
        cbegin = cumsum(x(1:width-2))
        cbegin = cbegin(1:2:end)./(1:2:(width-2))'
        cend = cumsum(x(n:-1:n-width+3))
        cend = cend(end:-2:1)./(width-2:-2:1)'
        y = [cbegin;y(width:end);cend];
        
    otherwise
        b = 'faulty method'
end


end

