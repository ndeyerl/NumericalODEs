% % Nicole Deyerl
% % MATH 6321 (Dan Reynolds)
% % 9/2/16
% % Homework 1, Problem 2 (Computational)
% % This script finds the roots of the lagrange interpolant of the points given
% % by problem 2 of homework 1.  It's a script which uses Newton's method to find
% % the roots of this particular function.

clear;
% lagrange interpolant of 4 data points
f = @(x) (-4/5 + 4/495 + 8/9)*x.^3 + (-28/5+8/495 + 52/9)*x.^2 +...
    (-31/5 + 13/1485 + 116/27)*x + (-119/135 + 14/13365 + 136/243); 
% derivative of lagrange interpolant
fprime = @(x) 3*(-4/5 + 4/495 + 8/9)*x.^2 + 2*(-28/5+8/495 + 52/9)*x +...
    (-31/5 + 13/1485 + 116/27);

% Newton's method (reference text: Numerical Analysis by Sauer from
% previous course)
n = 10; %number of steps to use in newton solver
x = zeros(1,n+1); %set up solution vector

% at most 3 roots -> 3 Newton solves with 3 initial guesses x0 = x(1)
x(1) = -1; % "x0"
for i = 1:n % start with "x0"
    x(i+1) = x(i) - f(x(i))/fprime(x(i)); % newton's method
end

fprintf('The first root is %d \n',x(n+1));

x(1) = -3;% "x0"
for i = 1:n
    x(i+1) = x(i) - f(x(i))/fprime(x(i)); 
end

fprintf('The second root is %d \n',x(n+1));

x(1) = -10;% "x0"
for i = 1:n
    x(i+1) = x(i) - f(x(i))/fprime(x(i));
end

fprintf('The third root is %d \n',x(n+1));