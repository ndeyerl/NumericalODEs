function [f] = fkaptrue(t)
%Analytic solution for Kaps problem
f = zeros(2,1);
f(1) = exp(-2*t);
f(2) = exp(-t);
end

