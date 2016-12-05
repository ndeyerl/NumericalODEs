function [f] = vdpfull(t,y,eps)
%Non-split VDP problem
%Note IC is [2,-0.6666654321121172] at t=0
f = [y(2); (1/eps)*((1-y(1).^2)*y(2) - y(1))];
end

%[t,y] = ode15s(@(t,y) vdpfull(t,y,eps),[0 10],[2,-0.6666654321121172]);
%plot(t,y(:,1))