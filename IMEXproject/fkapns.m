function [f] = fkapns(y,t)
%Nonstiff part of Kaps' problem
%Note IC is [1,1] at t=0
f = zeros(2,1);
f(1) = -2*y(1);
f(2) = y(1)-y(2)-y(2).^2;
end

