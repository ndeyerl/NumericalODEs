function [z] = ARKERKDIRKstages(eps, s, AE, AI, c, h, un, t, pnum, sz)
z = zeros(round(sz),round(s));
z(:,1) = un;
%guess = zeros(2,1);
if(pnum==1)%kaps problem
    nsfunc = @(z,t)fkapns(z,t);
    sfunc = @(e,z,t)fkaps(e,z,t);
    residfunc = @(e,h,A,u,sum,t)fkapresid(e,h,A,u,sum,t);
elseif(pnum==2) %prothero robinson trig problem
    nsfunc = @(z,t)fprtrigns(z,t);
    sfunc = @(e,z,t)fprtrigs(e,z,t);
    residfunc = @(e,h,A,u,sum,t)fprtrigresid(e,h,A,u,sum,t);
elseif(pnum==3) %hw3 prob 5 variable stiffness scalar ode
    nsfunc = @(z,t)fhwns(z,t);
    sfunc = @(eps,z,t)fhws(eps,z,t);
    residfunc = @(e,h,A,u,sum,t)fhwresid(e,h,A,u,sum,t);
elseif(pnum==4) %(8.8) scalar ode
    nsfunc = @(z,t)f88ns(z,t);
    sfunc = @(eps,z,t)f88s(eps,z,t);
    residfunc = @(e,h,A,u,sum,t)f88resid(e,h,A,u,sum,t);
elseif(pnum==200) %vdp problem
    nsfunc = @(z,t)fvdpns(z,t);
    sfunc = @(e,z,t)fvdps(e,z,t);
    residfunc = @(e,h,A,u,sum,t)fvdpresid(e,h,A,u,sum,t);
end
for i = 2:s
    %fprintf('stage loop\n');
    sum = zeros(round(sz),1);
    for j = 1:i-1
        fns = nsfunc(z(:,j),t+c(j)*h);
        fs = sfunc(eps,z(:,j),t+c(j)*h);
        sum = sum + AE(i,j)*fns + AI(i,j)*fs;
    end
    resid = residfunc(eps,h,AI(i,i),un,sum,t+c(i)*h);
    fns = nsfunc(z(:,i-1),t+c(i)*h);
    fs = sfunc(eps,z(:,i-1),t+c(i)*h);
    guess = z(:,i-1) + h*fns + h*fs;
    options = optimoptions('fsolve','Display','none');
    z(:,i) = fsolve(resid,guess,options);
%    z(:,i) = fsolve(resid,z(:,i-1));
end
end