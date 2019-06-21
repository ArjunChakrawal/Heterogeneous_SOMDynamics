function ks_mult=ks_mult_estimation(T,cs0,cb0,co20,kb,ksmm, km, Imic,Y,rho, Vol_mic)

options = odeset('Stats','off','AbsTol',1e-6,'RelTol',1e-6);
nstep=500; t=linspace(0,T,nstep);

St=ones(1, length(t)).*Imic;
func=@(tt,c)[-ksmm *c(1)*c(2)/(km+c(1)) + kb*c(2) + interp1(t, St,tt);...
    Y*ksmm *c(1)*c(2)/(km+c(1)) - kb*c(2); (1-Y)*ksmm*c(1)*c(2)/(km+c(1))];
[~, Ch]=ode45(func,t,[cs0,cb0,co20], options); %homogeneous solution 

ks=(ksmm.*Y - kb)./(Y.*km);
p0=ks;
optionsfmin = optimset( 'MaxFunEvals', 1000,'MaxIter', 1000, 'TolFun', 1e-10, 'TolX',1e-10);
p_estimate  = fminsearch(@(p)ode_fitting(t,Ch,p,St),p0,optionsfmin);
ks_mult=p_estimate./(0.001*rho.*1e-12.*1e15*Vol_mic);

    function err1 = ode_fitting(exp_t,exp_y,p,St)
        f=@(tt,c)[-p *c(1)*c(2) + kb*c(2) + interp1(exp_t, St,tt);...
                Y*p *c(1)*c(2) - kb*c(2); (1-Y)*p*c(1)*c(2)];
            [~, y]=ode45(f,exp_t,[exp_y(1,1),exp_y(1,2),exp_y(1,3)], options);
        err1= sqrt(mean((y(:,2)-exp_y(:,2)).^2));      
    end

end