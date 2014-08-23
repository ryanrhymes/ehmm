function [x,y] = HMMsample(pi, Q, g, n)
%
%  in  : pi = initial distribution of the HMM states
%        Q = transition matrix
%        n = number of samples
%  out : [x,y] = sample trajectory of size n of a HMM defined by (pi, Q, g):
%        x = sample trajectory of size n of a Markov Chain with initial distribution pi and transition matrix Q
%        y = observations such that the conditionnal distribution on x, i.e. given x(k) is g(x(k), :)
%

cQ = cumsum(Q, 2);
x = zeros(1, n); x(1) = 1+sum(rand>cumsum(pi));
y = zeros(1, n); y(1) = exprnd(g(x(1)));

for j=2:n
    x(j) = 1+sum(rand>cQ(x(j-1),:));
    y(j) = exprnd(g(x(j)));
end