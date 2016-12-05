function [solution] = Newton(MyFunc,Jacobian,Guess,tol)
% solves the non-linear vector equation F(x)=0
% set using Newton Raphson iteration
% INPUTS;
% Myfunc=Handle to the function which returns the vectorF(x)
% Jacobian=Handle to the function which returns the JacobianMatrix
% Guess = Initial Guess (a vector)
% tol = Tolerance
%
% OUTPUTS
% solution = The solution to F(x) = 0
% for tests:
% F = feval(MyFunc,g); error = max(abs(F));  J = feval(Jacobian,g);dx = J\(-F);    g = g+dx;
g = Guess;
%x = sym('x',[2 1]);
%set the error 2*tol to make sure the loop runs at least once 
error = 2*tol;
counter = 0;
while error > tol
    %disp(error)
    %calculate the function values at the current iteration
    F = feval(MyFunc,g);
    %F = subs(MyFunc,x(1),g(1));
    %F = subs(F,x(2),g(2));
    %F = double(F);
    %calculate the error
    error = max(abs(F));
    %calculate the jacobian matrix
    J = feval(Jacobian,g);
    [L,U] = ilu(sparse(J));
    %J = subs(Jacobian,x(1),g(1));
    %J = subs(J,x(2),g(2));
    %J = double(J);
    %calculate the update (solve the linear system)
    %dx = J\(-F);
    Y = L\(-F);
    dx = U\Y;

    %update the x value
    g = g+dx;
    counter = counter+1;
end %while loop
solution = g;

return