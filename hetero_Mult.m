function hetero_Mult(cs,cb,sav_fname,T,t,Imic,ks, kb, Y, nx, ny, options)

hbar = parfor_progressbar(100,'Parallel loop from hetero Mult is running, please wait!...');  %create the progress bar

parfor i=1:nx
    hbar.iterate(1);   % update progress by one iteration
    for j=1:ny
        St=ones(length(t), 1).*Imic;
        f=@(tt,c)[-ks(i,j) *c(1)*c(2) + kb(i,j)*c(2) + interp1(t, St,tt);...
            Y*ks(i,j) *c(1)*c(2) - kb(i,j)*c(2); (1-Y)*ks(i,j)*c(1)*c(2)];
        sol=ode45(f,[0,T],[cs(i,j),cb(i,j),0], options);
        C=(deval(sol, t))';
        Cs{i,j}=C(:,1);
        Cb{i,j}=C(:,2);
        CO2{i,j}=C(:,3);
    end       
end
delete(hbar)
% allocating differnt variables
Cs_mean=zeros(length(t),1);
Cs_std=Cs_mean;Cs_m3=Cs_mean;
Cb_std=Cs_mean;
CO2_mean=Cs_mean;CO2_std=Cs_mean;
dec_m=Cs_mean;dec_std=Cs_mean;
cov_cscb=Cs_mean;
cov_csks=cov_cscb;
cov_cbks=cov_cscb;
Cs_m4=zeros(length(t),1);
Cb_mean=zeros(length(t),1);
rho_CsCb=zeros(length(t),1);
rho_Csks= rho_CsCb;
rho_Cbks=rho_CsCb;
E_cs2_cb=zeros(length(t),1);
E_cs_cb2=zeros(length(t),1);
E_cs2_cb2=zeros(length(t),1);
E_cs_cb_ks=zeros(length(t),1);
s11=zeros(nx,ny);s22=s11;s33=s11;dec=s33;
% Figures to plot movies of cs and cb
% figure;
% hax1=axes;
% figure;
% hax2=axes;

for k=1:length(t)
 
    for i=1:nx
        for j=1:ny
            s1=Cs(i,j);
            s1=cell2mat(s1);
            s11(i,j)=s1(k);
            s2=Cb(i,j);
            s2=cell2mat(s2);
            s22(i,j)=s2(k);
            s3=CO2(i,j);
            s3=cell2mat(s3);
            s33(i,j)=s3(k);        
            dec(i,j)= ks(i,j)*s11(i,j)*s22(i,j);
        end
    end

CO2_mean(k)=mean2(s33);
CO2_std(k)=std2(s33);

dec_m(k)=mean2(dec);
dec_std(k)=std2(dec);
Cs_mean(k)=mean2(s11);
Cs_std(k)=std2(s11);
Cs_m4(k)=moment2(s11,4);
Cs_m3(k)=moment2(s11,3);

Cb_mean(k)=mean2(s22);
Cb_std(k)=std2(s22);

rho_CsCb(k)=corr2(s11,s22);
rho_Csks(k)=corr2(s11,ks);
rho_Cbks(k)=corr2(s22,ks);

cov_cscb(k)= rho_CsCb(k)*std2(s11)*std2(s22);
cov_csks(k)= rho_Csks(k)*std2(s11)*std2(ks);
cov_cbks(k)= rho_Cbks(k)*std2(s22)*std2(ks);

E_cs2_cb(k) = exp_A2_B(s11,s22);
E_cs_cb2(k) = exp_A_B2(s11,s22);
E_cs2_cb2(k) = exp_A2_B2(s11,s22);
E_cs_cb_ks(k) = exp_CsCbks(s11,s22,ks);

s11=zeros(nx,ny);
s22=zeros(nx,ny);
k
end
save(sav_fname)
end


