% disturbance testing


t = 0:0.01:10;
t0 = 5;
sig = 1
A = 1
gauss = @(t,sig,t0)  1/(sig*sqrt(2*pi)) *exp(-(t-t0).^2/sig^2/2) .* sin(pi/sig*t);


% sum(gaussian*0.01)
figure();
plot(t,gauss(t,sig,t0))
