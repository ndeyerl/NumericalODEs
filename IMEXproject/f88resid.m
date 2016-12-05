function [resid] = f88resid(eps,h,AI,un,sum,t)
%Residual for trigonometric prothero robinson problem stage calculations
resid = @(x) x-un-h*sum-h*AI*((1/eps)*(x-cos(t)-sin(t)));
end

