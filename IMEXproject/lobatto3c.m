function [u,t] = lobatto3c(A, b, c, h, range, eps, pnum)
if(pnum==200)
    myfunc = @(t,y,eps)vdpfull(t,y,eps);
    ic = [2,-0.6666654321121172];
end
n = (range(2)-range(1))/h; %number of steps
u = zeros(2,n); %2 rows by n colns->holds solution
t = zeros(1,n);
u(:,1) = ic; %set initial condition
t(1) = range(1); %set initial time
s = 3; %3 stages
for k=1:n-1
    %fprintf('main loop\n');
    t(k+1) = t(k) + h;
    %fprintf('t = %d\n',t);
    z = lobatto3c(eps, A, c, h, u(:,k), t(k+1), pnum);
    sum = zeros(2,1);
    for i=1:s
        f = myfunc(t(k+1),z(:,i),eps);
        sum = sum + b(i)*f;
    end
    u(:,k+1) = u(:,k) + h*sum;
end

end