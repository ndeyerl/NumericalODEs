function [f] = fhws(eps,y,t)
%Stiff part of hw3 prob 5 variable stiffness problem
%Note IC is [0] at t=0
f = (1/eps)*y -(1/eps)*atan(t);
end

