function [resid] = fhwresid(eps,h,AI,un,sum,t)
%Residual for hw3 prob 5 variable stiffness problem stage calculations
resid = @(x) x-un-h*sum-h*AI*((1/eps)*x-(1/eps)*atan(t));
end

