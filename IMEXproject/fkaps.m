function [f] = fkaps(eps,y,t)
%Stiff part of Kaps' problem
%Note IC is [1,1] at t=0
f = zeros(2,1);
f(1) = -(1/eps)*y(1) + (1/eps)*y(2).^2;
f(2) = 0.0;
end

