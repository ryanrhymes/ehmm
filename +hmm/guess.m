function [ pi, Qg, gg ] = guess( states )
% Return the initial guess on HMM model parameters

pi = rand(1,states); pi = pi / sum(pi);

if (states==2)
    Qg = [.5 .5; .5 .5];
    gg = [1; 0.5];
elseif (states==3)
    Qg = [.3 .3 .4; .3 .4 .3; .4 .3 .3];
    gg = [1.6; .4; .1];
elseif (states==4)
    Qg = repmat([.25 .25 .25 .25],4,1);
    gg = [6.4; 1.6; .4; .1];
end

end

