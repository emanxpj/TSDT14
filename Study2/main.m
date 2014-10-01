addpath ./Study1;

N = 2^10;
N2 = N+1;
Ts = 1; %length of the measured signal.
fs = N/Ts; %sampling frequency.
T = Ts/N; %sampling length.

x = randn(1,N);


n = linspace(0,N,N); 
wc = 0.4;
[b,a] = butter(10,wc,'low');

y = filter(b,a,x);