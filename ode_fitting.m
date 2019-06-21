function err1 = ode_fitting(exp_t,exp_y,p,kb,St)
options = odeset('Stats','off','AbsTol',1e-6,'RelTol',1e-6);
Y=0.31;
[~,y] = DEC_multiplicative(exp_y(1,1),exp_y(1,2),exp_y(1,3),exp_t, p(1), kb, Y,St, options);
err1= sqrt(mean((y(:,3)-exp_y(:,3)).^2));
end
