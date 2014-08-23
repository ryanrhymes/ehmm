function [Q, g, l] = HMMbaumwelch(y, pi, tol, maxIt, Q, g)
%
%  in : y = observation vector
%       pi = initial distribution of the states
%       tol = tolerance for the stopping criterion
%       maxIt = maximal number of iterations
% out : Q = estimate of the transition matrix
%       g = estimated probabilities of transition: gh(x,y) = estimate of P(Y=y | X=x) for 1<=x<=k
%       l = log-likelihood of y for parameters Q and g
%

import hmm.*;

k = length(pi); n = length(y);
it = 0; oldQ = Q; oldg = g+tol+1;
while ((norm(oldQ(:)-Q(:), 1) + norm(oldg-g, 1) > tol) && (it<maxIt))
  disp([it, norm(oldQ(:)-Q(:), 1) + norm(oldg-g, 1)]);

  it = it + 1;
  % calculate the posterior given the model
  [phi, c] = HMMfilter(y, pi, Q, g);
  beta = HMMsmoother(y, Q, g, c);
  post = phi .* beta;
  
  % e-step: average transitions given the model
  h = cell2mat(arrayfun(@(x) exppdf(x, g), y(2:end), 'UniformOutput', false));
  N =Q.*(phi(:, 1:(end-1))*(beta(:, 2:end).*h./(ones(k, 1)*c(2:end)))');   
  % e-step: average durations
  yt = repmat(y, k, 1);
  Ma = sum(post .* yt, 2); Mb = sum(post, 2);
  
  % m-step: maximize the log-likelihood and update model parameters
  oldQ = Q; oldg = g;
  Q = N ./ (sum(N, 2) * ones(1, k));
  g = Ma ./ Mb;
end
l = sum(log(c));