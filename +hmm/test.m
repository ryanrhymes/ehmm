% TEST Baum-Welch Algorithm
clc; close all; clear all;
import hmm.*;

n = 1000;
% 2-state HMM
pi = [.5 .5];
Q = [0.9, 0.1; 0.3, 0.7];
g = [1.85; 0.25];
Qg = [.5 .5; .5 .5];
gg = [1; 0.5];
% 3-state HMM
% pi = [.3 .3 .4];
% Q = [.8 .1 .1; .3 .3 .4; .5 .1 .4];
% g = [1.85; 0.55; 0.20];
% Qg = [.4 .3 .3; .3 .3 .4; .6 .3 .1];
% gg = [1; .8; .5];

[x,y] = HMMsample(pi, Q, g, n);
[phi, c] = HMMfilter(y, pi, Q, g);
[~, xh] = max(phi, [], 1);
mean(x==xh)  % estimation performance of the filter

beta = HMMsmoother(y, Q, g, c);
psi = phi .* beta; % posterior probabilities
[~, xh] = max(psi, [], 1); % MAP
mean(x==xh) % estimation performance of the marginal MAP

[Qh, gh, lh] = HMMbaumwelch(y, pi, 1e-16, 1000, Qg, gg);
Q, Qh, g, gh

% compare real and maximum log-likelihood
[phi, ch] = HMMfilter(y, pi, Qh, gh);
sum(log(ch))-sum(log(c))
[~, xh] = max(phi, [], 1);
mean(x==xh)