function [f] = fvdpns(y,t)
%Nonstiff part of VDP problem
%Note IC is [2,?0.6666654321121172] at t=0
f = zeros(2,1);
f(1) = y(2);
f(2) = 0.0;
end

