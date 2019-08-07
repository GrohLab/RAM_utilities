

function S=makeFilteredStim(N,dt,tau,sigma,mu,rseed)
%function S=makestim(N,dt,tau,sigma)
%%produces Gaussian stimulus N samples long, with std = sigma
%%1/sample rate =dt
%%rseed resets random number generation
randn('state',rseed);
I = sigma*randn(N,1)+mu; % stimulus
S = zeros(N,1);
for i = 1:N
   S(i+1) = (S(i) + I(i))/(1+dt/tau);  % exponential filtering
end
S=S(2:end);
S=S./std(S)*sigma;
