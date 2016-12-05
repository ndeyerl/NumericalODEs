function [f] = fvdps(eps,y,t)
%Stiff part of VDP problem
%Note IC is [2,-0.6666654321121172] at t=0
f = zeros(2,1);
f(1) = 0.0;
f(2) = (1/eps)*((1-y(1).^2)*y(2) - y(1));
end

