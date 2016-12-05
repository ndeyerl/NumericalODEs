function [z] = lobatto3cstages(eps, A, c, h, un, t, pnum)
z = zeros(2,s);
z(:,1) = un;
%guess = zeros(2,1);
if(pnum==200) %vdp problem
    myfunc = @(t,y,eps)vdpfull(t,y,eps);
    residfunc = @(e,h,A,u,sum,t)fvdpresid(e,h,A,u,sum,t);
end
for i = 2:s
    %fprintf('stage loop\n');
    sum = zeros(2,1);
    for j = 1:i-1
        f = myfunc(eps,z(:,j),t);
        sum = sum + AE(i,j)*fns + ...
            AI(i,j)*fs;
    end

%     resid = sym([x(1)-un(1)-h*sum(1)-h*AI(i,i)*((-1/eps)*x(1)+(1/eps)*x(2).^2);...
%                        x(2)-un(2)-h*sum(2)-h*AI(i,i)*0]);
%     residJac = sym([1-h*AI(i,i)*(-1/eps), 1-h*AI(i,i)*2*(1/eps)*x(2); 1,1]);
    resid = residfunc(eps,h,AI(i,i),un,sum,t);
%    resid = @(x) [x(1)-un(1)-h*sum(1)-h*AI(i,i)*((-1/eps)*x(1)+(1/eps)*x(2).^2);...
%                        x(2)-un(2)-h*sum(2)-h*AI(i,i)*0];
    %residJac = @(x) [1-(-h*AI(i,i)*(-1/eps)), -(h*AI(i,i)*2*(1/eps)*x(2));...
    %                 0,1];
    %z(:,i) = Newton(resid,residJac,z(:,i-1),10.^(-11),25);
    %z(:,i) = fsolve(resid,z(:,i-1));
    %fprintf('start newton\n');
    %z(:,i) = newtonm(z(:,i-1),resid,residJac);
    fns = nsfunc(z(:,i-1),t);
    fs = sfunc(eps,z(:,i-1),t);
    guess = z(:,i-1) + h*fns + h*fs;
%    z(:,i) = Newton(resid,residJac,z(:,i-1),10.^-12);
%    options = optimoptions('fsolve','Display','iter');
    z(:,i) = fsolve(resid,guess);
    
    %z(:,i) = Newton(resid,residJac,z(:,i-1),10.^-5);
    %fprintf('newton completed\n');
end
end