function [f] = fprtrigs(eps,y,t)
%Stiff part of trigonometric prothero robinson problem
%Note IC is [1] at t=0
f = (1/eps)*(y-cos(t));
end

