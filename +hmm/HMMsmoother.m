function beta = HMMsmoother(y, Q, g, c)
%
%  in : y = observation vector (1xn), 1<=y<=r
%       Q = transition matrix of (kxk)
%       g = emission matrix (kx1)
%       c = vector computed by function filter: c(t) = P(Y(t) = y(t)| Y(1:t-1)=y(1:t-1))
% out : beta = smoothing factors: 
%       beta(x,t) = P(Y(t+1:n)=y(t+1:n) | X(t)=x) / P(Y(t+1:n)=y(t+1:n) | Y(1:t)=y(1:t)) for 1<=x<=k and 1<=t<=n
%       post = phi .* beta; i.e. post(x,t) = P(X(t)=x | Y(1:n)=y(1:n))
%

n = length(y);
beta = ones(size(Q, 1), n);
for t=(n-1):-1:1
  beta(:, t) = Q * (exppdf(y(t+1),g) .* beta(:,t+1)) / c(t+1);
end