function [phi, c] = HMMfilter(y, pi, Q, g)
%
%  in : y = vector of observations of size n, with values between 1 and r
%       pi = initial distribution as vector of size k
%       Q = transition matrix of size k
%       g = emission matrix of size k x r
% out : phi = filter(x,t) = P(X(t)=x | Y(1:t)=y(1:t)) for 1<=x<=k and 1<=t<=n
%       c(t) = P(Y(t) = y(t)| Y(1:t-1)=y(1:t-1))
%


n = length(y); k = length(g);
phi = zeros(size(Q, 1), n);
c = zeros(1, n);

alpha = pi'.*exppdf(y(1),g);
c(1)= sum(alpha);
phi(:, 1) = alpha/c(1);

for t=2:n
  alpha = (phi(:, t-1)' * Q)' .* exppdf(y(t),g);
  c(t) = sum(alpha);
  if c(t) > 0
    phi(:, t) = alpha / c(t);
  else
    c(t) = rand();
    vec = rand(k, 1);
    phi(:, t) = vec / sum(vec);
  end
end