function [t,w]=RK4(inter,y0,n,eps)
t(1)=inter(1);
w(:,1)=y0;
h=(inter(2)-inter(1))/n;
for i=1:n;
    t(i+1)=t(i)+h;
    w(:,i+1)=RK4step(t(i),w(:,i),h,eps);
end
end

function w=RK4step(t,w,h,eps)
k1=yprime(t,w,eps);
k2=yprime(t+h/2,w+(h/2)*k1,eps);
k3=yprime(t+h/2,w+(h/2)*k2,eps);
k4=yprime(t+h,w+h*k3,eps);
w=w+(h/6)*(k1+2*k2+2*k3+k4);
end



% 
% h = (b-a)/N;        %the step size
% t(1) = a;
% w(:,1) = alpha;     %initial conditions
% 
% for i = 1:N
%    k1 = h*f(t(i), w(:,i));
%    k2 = h*f(t(i)+h/2, w(:,i)+0.5*k1);
%    k3 = h*f(t(i)+h/2, w(:,i)+0.5*k2); 
%    k4 = h*f(t(i)+h, w(:,i)+k3);
%    w(:,i+1) = w(:,i) + (k1 + 2*k2 + 2*k3 + k4)/6;
%    t(i+1) = a + i*h;
% end
% 
% 



function z=yprime(t,w,eps)
z(1)=w(2); 
z(2)=(1/eps)*((1-w(1).^2)*w(2) - w(1));
end