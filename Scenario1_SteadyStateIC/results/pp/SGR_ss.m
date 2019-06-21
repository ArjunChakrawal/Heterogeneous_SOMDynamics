clc
close all
clearvars;

ksmm=0.018;
kb=0.00028;  %h-1
Y=0.31;
km=25;
rho= 1.65; %g/cm3 soil bulk density 
rho_wood=1.1;%g/cm3  max carboff density in form of wood 
nx=100;ny=100; %offe pixel = 25microff

Vol_mic = 50*50*50; %um^3
Vmic_cc=Vol_mic*1e-12; %cm3
Vdomain_cc=Vol_mic*nx*ny*1e-12 ; % cm^3

selpath = uigetdir;
addpath(selpath)
addpath([selpath,'\Scenario1_SteadyStateIC\results\'])
LC=linspecer(3); % for disctinctive line color 

%% homo
load('Ph_mult_neg_ss_1.mat')
cs0=mean2(cs);cb0=mean2(cb);co20=0;
ks=mean2(ks);
kb=mean2(kb);
c(1,:)=[cs0,cb0,co20];
St=ones(1, length(t)).*Imic;
[~,C1]= DEC_multiplicative(cs0,cb0,co20,t, ks, kb, Y,St, options);
dec0=ks.*C1(:,1).*C1(:,2);
%% Figure properties
phi = (1+sqrt(5))/2; %godlen ratio
width = 6;     % Width in inches
height = width/phi;    % Height in inches
alw = 0.5;    % AxesLineWidth
fsz = 10;      % FontSize
fsz_title= fsz+1;
fsz_ax_lgnd = fsz-1;
lw = 1.0;      % LineWidth
msz = 8;       % MarkerSize
AHlen=5;
AHwid=4;
AHlw=0.01;
hstyle='vback3';
clr = 'k';     % Line color
ftype1= 'Helvetica';  %hLegend, gca
ftype2='Times New Roman'; % hTitle, hXLabel, hYLabel
xycolor= [0.2 0.2 0.2];
fig=figure;
pos = get(gcf, 'Position');
set(gcf, 'Position', [1228 267 width*100, height*100]); %<- Set size
% [ha, pos] = tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right])
[ha, pos] = tight_subplot(1, 2, [.035 .05],[.15 .15],[.15 .05]);
a1=ha(1);a2=ha(2);

TT=10;
xtick= 5:0.5:7.5;
ytick= 5.25:0.5:7.25;
aa=25;
% cstck=[0,40,80];
%%

load('Ph_mult_pos_ss_1.mat');
csm = Cs_mean.*1e-12./(rho*Vmic_cc);
sgr=((1-Y).*dec_m./Cb_mean);
plot(a1,csm, sgr,'--','linewidth',lw, 'Color', clr); hold(a1,'on'); 
xlim(a1,[xtick(1), xtick(end)])
ylim(a1,[ytick(1), ytick(end)].*1e-4)
csmt=csm(1:aa);
sgrt=sgr(1:aa);
axes(a1);
for i=1:length(csmt)-1
x_begin=csmt(i);
x_end=csmt(i+1);
y_begin=sgrt(i);
y_end= sgrt(i+1);
xa= [x_begin x_end];ya= [y_begin y_end];
[xaf,yaf] = ds2nfu(xa,ya); % Convert to normalized figure units 
A=annotation('arrow',xaf,yaf); % Add annotation
A.HeadLength=AHlen;
A.HeadWidth=AHwid;
A.LineStyle='--';
A.LineWidth=AHlw;
A.HeadStyle=hstyle;
A.Color=LC(1,:);
end 


load('Ph_mult_neg_ss_1.mat');
csm = Cs_mean.*1e-12./(rho*Vmic_cc);
sgr=((1-Y).*dec_m./Cb_mean);
plot(a1,csm, sgr,'-.','linewidth',lw, 'Color', clr);  
xlim(a1,[xtick(1), xtick(end)])
ylim(a1,[ytick(1), ytick(end)].*1e-4)

csmt=csm(1:aa);
sgrt=sgr(1:aa);
axes(a1);
for i=1:length(csmt)-1
x_begin=csmt(i);
x_end=csmt(i+1);
y_begin=sgrt(i);
y_end= sgrt(i+1);
xa= [x_begin x_end];ya= [y_begin y_end];
[xaf,yaf] = ds2nfu(xa,ya); % Convert to normalized figure units 
A=annotation('arrow',xaf,yaf); % Add annotation
A.HeadLength=AHlen;
A.HeadWidth=AHwid;
A.LineStyle='none';
A.LineWidth=AHlw;
A.HeadStyle=hstyle;
A.Color=LC(2,:);

end 

load('Ph_mult_zero_ss_1.mat');
csm = Cs_mean.*1e-12./(rho*Vmic_cc);
sgr=((1-Y).*dec_m./Cb_mean);
plot(a1,csm, sgr,':','linewidth',lw, 'Color', clr); 
xlim(a1,[xtick(1), xtick(end)])
ylim(a1,[ytick(1), ytick(end)].*1e-4)

csmt=csm(1:aa);
sgrt=sgr(1:aa);
axes(a1);
for i=1:length(csmt)-1
x_begin=csmt(i);
x_end=csmt(i+1);
y_begin=sgrt(i);
y_end= sgrt(i+1);
xa= [x_begin x_end];ya= [y_begin y_end];
[xaf,yaf] = ds2nfu(xa,ya); % Convert to normalized figure units 
A=annotation('arrow',xaf,yaf); % Add annotation
A.HeadLength=AHlen;
A.HeadWidth=AHwid;
A.LineStyle='none';
A.LineWidth=AHlw;
A.HeadStyle=hstyle;
A.Color=LC(3,:);

end

scatter(a1,C1(1,1).*1e-12./(rho*Vmic_cc), ((1-Y).*dec0(1)./C1(1,2)),60,'filled','ko','linewidth',lw);

%%

load('Ph1ch2_mult_pos_ss.mat')
csm = Cs_mean.*1e-12./(rho*Vmic_cc);
sgr=((1-Y).*dec_m./Cb_mean);
plot(a2,csm, sgr,'--','linewidth',lw, 'Color', clr); hold(a2,'on'); 
xlim(a2,[xtick(1), xtick(end)])
ylim(a2,[ytick(1), ytick(end)].*1e-4)

csmt=csm(1:aa);
sgrt=sgr(1:aa);
axes(a2);
for i=1:length(csmt)-1
x_begin=csmt(i);
x_end=csmt(i+1);
y_begin=sgrt(i);
y_end= sgrt(i+1);
xa= [x_begin x_end];ya= [y_begin y_end];
[xaf,yaf] = ds2nfu(xa,ya); % Convert to normalized figure units 
A=annotation('arrow',xaf,yaf); % Add annotation
A.HeadLength=AHlen;
A.HeadWidth=AHwid;
A.LineStyle='none';
A.LineWidth=AHlw;
A.HeadStyle=hstyle;
A.Color=LC(1,:);

end

load('Ph1ch2_mult_neg_ss.mat')
csm = Cs_mean.*1e-12./(rho*Vmic_cc);
sgr=((1-Y).*dec_m./Cb_mean);
plot(a2,csm, sgr,'-.','linewidth',lw, 'Color', clr);  
xlim(a2,[xtick(1), xtick(end)])
ylim(a2,[ytick(1), ytick(end)].*1e-4)

csmt=csm(1:aa);
sgrt=sgr(1:aa);
axes(a2);
for i=1:length(csmt)-1
x_begin=csmt(i);
x_end=csmt(i+1);
y_begin=sgrt(i);
y_end= sgrt(i+1);
xa= [x_begin x_end];ya= [y_begin y_end];
[xaf,yaf] = ds2nfu(xa,ya); % Convert to normalized figure units 
A=annotation('arrow',xaf,yaf); % Add annotation
A.HeadWidth=AHwid;
A.HeadLength=AHlen;
A.LineStyle='none';
A.LineWidth=AHlw;
A.HeadStyle=hstyle;
A.Color=LC(2,:);

end

load('Ph1ch2_mult_zero_ss.mat')
csm = Cs_mean.*1e-12./(rho*Vmic_cc);
sgr=((1-Y).*dec_m./Cb_mean);
plot(a2,csm, sgr,':','linewidth',lw, 'Color', clr);  
xlim(a2,[xtick(1), xtick(end)])
ylim(a2,[ytick(1), ytick(end)].*1e-4)

csmt=csm(1:aa);
sgrt=sgr(1:aa);
axes(a2);
for i=1:length(csmt)-1
x_begin=csmt(i);
x_end=csmt(i+1);
y_begin=sgrt(i);
y_end= sgrt(i+1);
xa= [x_begin x_end];ya= [y_begin y_end];
[xaf,yaf] = ds2nfu(xa,ya); % Convert to normalized figure units 
A=annotation('arrow',xaf,yaf); % Add annotation
A.HeadLength=AHlen;
A.HeadWidth=AHwid;
A.LineStyle='none';
A.LineWidth=AHlw;
A.HeadStyle=hstyle;
A.Color=LC(3,:);

end

S1=scatter(a2,C1(1,1).*1e-12./(rho*Vmic_cc), ((1-Y).*dec0(1)./C1(1,2)),60,'filled','ko','linewidth',lw);
xlim(a2,[xtick(1), xtick(end)])

%%
ylabel(a1,'$SGR$ \ $\mathrm{(h^{-1})}$','Interpreter','latex','FontSize', fsz, 'FontName', ftype2);

hxlable=xlabel(a1,'$\overline{C}_S$ (mgC/gSoil)','Interpreter','latex','FontSize', fsz, 'FontName', ftype2);
hx=hxlable.Position;
hxlable.Position= [ hx(1)+0.85, hx(2)-1e-5,hx(3)];


t1=title(a1, 'Biophysical heterogeneity','FontSize', fsz_title, 'FontName', ftype2,'FontWeight' , 'normal', 'Interpreter','latex' );
titlePos1 = get( t1 , 'position');
titlePos1(2) = 7.4e-4;
set( t1 , 'position' , titlePos1);

t2=title(a2, 'Full heterogeneity','FontSize', fsz_title, 'FontName', ftype2,'FontWeight' , 'normal', 'Interpreter','latex', 'Interpreter','latex' );
titlePos2 = get( t2 , 'position');
titlePos2(2) = 7.4e-4;
set( t2 , 'position' , titlePos2);

hlegend=legend(a1,'Positive', 'Negative','Uncorrelated','Homogeneous');
set(hlegend,'box', 'off','Location', 'southeast', 'Interpreter','latex','FontSize', fsz_ax_lgnd, 'FontName', ftype1);

set(a1,'XTick',xtick,'YTick',ytick.*1e-4,'FontSize', fsz,'FontName','Times New Roman', 'LineWidth', alw, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on','XColor', xycolor, 'YColor', xycolor);
set(a2,'XTick',xtick,'YTick',ytick.*1e-4,'YTickLabel',[],'FontSize', fsz,'FontName','Times New Roman', 'LineWidth', alw, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on','XColor', xycolor, 'YColor', xycolor);
ytickformat(a1,'%1.2f')
str={'a','b'};
for i=1:2
    axes(ha(i)); %set the current axes to axes2
    text(ha(i).XLim(2)-ha(i).XLim(2)/30,  ha(i).YLim(2) - ha(i).YLim(2)/50, ['(',str{i},')'], 'FontWeight','normal','FontSize', fsz+2, 'FontName', ftype2)
end

% save figure 
set(gcf, 'Color', 'w')
export_fig(gcf,'SGR_ss.pdf','-r300')

%%
rmpath([selpath,'\Scenario1_SteadyStateIC\results\'])