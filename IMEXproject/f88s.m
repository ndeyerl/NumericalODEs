function [f] = f88s(eps,y,t)
%Stiff part of (8.8) problem
%Note IC is [1] at t=0
f = (1/eps)*(y-cos(t)-sin(t));
end

