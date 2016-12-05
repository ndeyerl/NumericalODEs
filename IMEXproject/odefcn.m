function dydt = odefcn(t,y,A,B)
dydt = zeros(2,1);
dydt(1) = y(2);
dydt(2) = (1/A)*(1-y(1).^2)*y(2) - y(1);
end
