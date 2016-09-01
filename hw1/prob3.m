% % Nicole Deyerl
% % MATH 6321 (Dan Reynolds)
% % 9/2/16
% % Homework 1, Problem 3
% % This script finds the root of a nonlinear system of equations given by 
% % homework 1.  It uses a relative tolerance of 10^06 and an absolute tolerance
% % of 10^-10, with an initial guess of (1,2).
% % % test change
%functions f1 and f2
f1 = @(x,y) x.^2 + y.^2 -4;
f2 = @(x,y) x*y - 1;

%Jacobian functions
Df1 = @(x,y) 2*x;
Df2 = @(x,y) 2*y;
Df3 = @(x,y) y;
Df4 = @(x,y) x;

%absolute tolerance threshold solver
abstol = 1; %initialize absolute tolerance
i = 1;
x = zeros(2,10); %solution matrix (made assuming while loop doesnt take more
                  %than 10 iterations
p = zeros(2,1); %intermediate solution vector
f = zeros(2,1); %function value vector
Df = zeros(2,2); %jacobian value vector

x(:,1) = [1,2]; %initial condition "x0"=(1,2)
while (abstol > 10.^(-10))    
    f = [ f1(x(1,i),x(2,i)); f2(x(1,i),x(2,i))]; %vector f(x(i))
    Df = [Df1(x(1,i),x(2,i)), Df2(x(1,i),x(2,i)); Df3(x(1,i),x(2,i)), ...
        Df4(x(1,i),x(2,i))]; %matrix Df(x(i))
    p(:,1) = Df\(-f); %solve Df*p =-f
    x(:,i+1) = x(:,i) + p(:,1); %solve for next x-value
    abstol = abs( (x(1,i+1)-x(1,i)) + (x(2,i+1)-x(2,i)) );
    reltol = abstol/abs(x(1,i) + x(2,i));
    i = i + 1;
end
% print solution
fprintf('--------------------\n');
fprintf('Absolute tolerance threshold = 10^-10\n');
fprintf('final abstol = %.3e \nsoln = (x,y) = (%.13d, %.13d)\n',...
    abstol, x(1,i),x(2,i)); %x0 has index 1, + i-1 time steps -> soln has index i


%relative tolerance threshold solver
reltol = 1; %initialize relative threshold
i = 1;
x = zeros(2,10); %re-initialize x
x(:,1) = [1,2]; %initial condition "x0"=(1,2)
while (reltol > 10.^(-6))  
    f = [ f1(x(1,i),x(2,i)); f2(x(1,i),x(2,i))]; %vector f(x(i))
    Df = [Df1(x(1,i),x(2,i)), Df2(x(1,i),x(2,i)); Df3(x(1,i),x(2,i)), ...
        Df4(x(1,i),x(2,i))]; %matrix Df(x(i))
    p(:,1) = Df\(-f); %solve Df*p =-f
    x(:,i+1) = x(:,i) + p(:,1); %solve for next x-value
    abstol = abs( (x(1,i+1)-x(1,i)) + (x(2,i+1)-x(2,i)) );
    reltol = abstol/abs(x(1,i) + x(2,i));
    i = i + 1;
end
%print solution
fprintf('--------------------\n');
fprintf('Relative tolerance threshold = 10^-6\n');
fprintf('final reltol = %.3e \nsoln = (x,y) = (%.13d, %.13d)\n',...
    reltol,x(1,i),x(2,i));
