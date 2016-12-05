function [resid] = fkapresid(eps,h,AI,un,sum,t)
%Residual for Kaps' problem stage calculations
x = sym('x',[2,1]);
resid = @(x) [x(1)-un(1)-h*sum(1)-h*AI*(-(1/eps)*x(1) + (1/eps)*x(2).^2);...
            x(2)-un(2)-h*sum(2)-h*AI*0.0];
end


