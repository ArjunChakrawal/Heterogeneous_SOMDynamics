function hetero_MM(cs,cb,sav_fname,T,t,Imic, ksmm, kb,km, Y, nx, ny, options)

d2dec_dcs2= @(csm,cbm,ksmmf,kmf)(-(2.*cbm.*kmf.*ksmmf)./(csm+ kmf).^3);
d2dec_dkm2= @(csm,cbm,ksmmf,kmf) (2*cbm*csm*ksmmf)./(csm + kmf).^3;
d2dec_dkmdksmm= @(csm,cbm,ksmmf,kmf) -(cbm*csm)./(csm + kmf).^2;
d2dec_dcsdcb= @(csm,cbm,ksmmf,kmf)(ksmmf.*kmf)./((csm+kmf).^2);
d2dec_dcsdksmm= @(csm,cbm,ksmmf,kmf)(cbm*kmf)./(csm + kmf).^2;
d2dec_dcbdksmm= @(csm,cbm,ksmmf,kmf)csm./(csm + kmf);
d2dec_dcsdkm= @(csm,cbm,ksmmf,kmf) (cbm*ksmmf*(csm - kmf))./(csm + kmf).^3;
d2dec_dcbdkm= @(csm,cbm,ksmmf,kmf) -(csm*ksmmf)./(csm + kmf).^2;

hbar = parfor_progressbar(100,'Parallel loop from hetero MM is running, please wait!');  %create the progress bar
parfor i=1:nx
    hbar.iterate(1);   % update progress by one iteration
        for j=1:ny
            St=ones(length(t), 1).*Imic;
            f=@(tt,c)[-ksmm(i,j) *c(1)*c(2)/(km(i,j)+c(1)) + kb(i,j)*c(2) + interp1(t, St,tt);...
                Y*ksmm(i,j) *c(1)*c(2)/(km(i,j)+c(1)) - kb(i,j)*c(2); (1-Y)*ksmm(i,j)*c(1)*c(2)/(km(i,j)+c(1))];
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
Cs_m4=zeros(length(t),1);
Cb_mean=zeros(length(t),1);
rho_CsCb=zeros(length(t),1);
E_cs2_cb=zeros(length(t),1);
E_cs_cb2=zeros(length(t),1);
E_cs2_cb2=zeros(length(t),1);
s11=zeros(nx,ny);s22=s11;s33=s11;dec=s33;

rho_cs_ksmm=rho_CsCb;
rho_cs_km=rho_CsCb;
rho_cb_ksmm=rho_CsCb;
rho_cb_km= rho_CsCb;
rho_ksmm_km=rho_CsCb;

cov_cs_ksmm=rho_cs_ksmm;
cov_cs_km=rho_CsCb;
cov_cb_ksmm=rho_CsCb;
cov_cb_km= rho_CsCb;
cov_ksmm_km=rho_CsCb;

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
            dec(i,j)= ksmm(i,j)*s11(i,j)*s22(i,j)/(km(i,j)+s11(i,j));
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
rho_cs_ksmm(k)=corr2(s11,ksmm);
rho_cs_km(k)=corr2(s11,km);
rho_cb_ksmm(k)=corr2(s22,ksmm);
rho_cb_km(k)= corr2(s22,km);
rho_ksmm_km(k)=corr2(ksmm,km);


cov_cscb(k)= rho_CsCb(k)*std2(s11)*std2(s22);
cov_cs_ksmm(k)= rho_cs_ksmm(k)*std2(s11)*std2(ksmm);
cov_cs_km(k)= rho_cs_km(k)*std2(s11)*std2(km);
cov_cb_ksmm(k)= rho_cb_ksmm(k)*std2(s22)*std2(ksmm);
cov_cb_km(k)= rho_cb_km(k)*std2(s22)*std2(km);
cov_ksmm_km(k)= rho_ksmm_km(k)*std2(ksmm)*std2(km);


cov1(k)=0.5*d2dec_dcs2(Cs_mean(k),Cb_mean(k),mean2(ksmm),mean2(km))*(std2(s22))^2;
cov2(k)=0.5*d2dec_dkm2(Cs_mean(k),Cb_mean(k),mean2(ksmm),mean2(km))*(std2(km))^2;
cov3(k)=d2dec_dkmdksmm(Cs_mean(k),Cb_mean(k),mean2(ksmm),mean2(km))*cov_ksmm_km(k);

cov4(k)=d2dec_dcsdcb(Cs_mean(k),Cb_mean(k),mean2(ksmm),mean2(km))*cov_cscb(k);
cov5(k)=d2dec_dcsdksmm(Cs_mean(k),Cb_mean(k),mean2(ksmm),mean2(km))*cov_cs_ksmm(k);
cov6(k)=d2dec_dcbdksmm(Cs_mean(k),Cb_mean(k),mean2(ksmm),mean2(km))*cov_cb_ksmm(k);
cov7(k)=d2dec_dcsdkm(Cs_mean(k),Cb_mean(k),mean2(ksmm),mean2(km))*cov_cs_km(k);
cov8(k)=d2dec_dcbdkm(Cs_mean(k),Cb_mean(k),mean2(ksmm),mean2(km))*cov_cb_km(k);


total_cov(k)= (1-Y).*(cov1(k) + cov2(k)+ cov3(k)+ cov4(k)+....
            cov5(k)+ cov6(k)+ cov7(k)+ cov8(k));


E_cs2_cb(k) = exp_A2_B(s11,s22);
E_cs_cb2(k) = exp_A_B2(s11,s22);
E_cs2_cb2(k) = exp_A2_B2(s11,s22);

k
s11=zeros(nx,ny);
s22=zeros(nx,ny);
 
end
close all
save(sav_fname)

