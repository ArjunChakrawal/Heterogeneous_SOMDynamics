function [t1,C]= DEC_multiplicative(cs0,cb0,co20, t, ks, kb, Y,St, options)
c(1,:)=[cs0,cb0,co20];
f=@(tt,c)[-ks *c(1)*c(2) + kb*c(2) + interp1(t, St,tt);...
    Y*ks *c(1)*c(2) - kb*c(2); (1-Y)*ks*c(1)*c(2)];
[t1, C]=ode45(f,t,c(1,:), options);
