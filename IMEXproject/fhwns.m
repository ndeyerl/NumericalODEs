function [f] = fhwns(y,t)
%Nonstiff part of hw3 prob 5 variable stiffness problem
%Note IC is [0] at t=0
f = 1/(1+t.^2);
end

