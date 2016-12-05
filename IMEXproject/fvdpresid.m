function [resid] = fvdpresid(eps,h,AI,un,sum,t)
%Residual for VDP problem stage calculations
x = sym('x',[2,1]);
resid = @(x) [x(1)-un(1)-h*sum(1)-h*AI*0.0;...
            x(2)-un(2)-h*sum(2)-h*AI*((1/eps)*((1-x(1).^2)*x(2) - x(1)))];
end

