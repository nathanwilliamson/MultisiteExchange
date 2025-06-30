%% Fig 11 source code
%% Williamson, et al, Magnetic Resonance Letters (2025)
%%%  figure showing representative raw exchange data from individual simulations.

%% simulate data
clear all
close all

set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

COLORSpatch= 1/255*[225,248,220; ... %light green. normal
                 254, 248, 221 ;...    % light yellow
                 255, 231, 199;...     % light orange. ouab, 10
                 247, 216, 186;...       % light red         
                 172,221,222;...   %greenish blue. sucrose
                 220, 225, 248;...  %purple. sucrose+ouab  
                 172-10,221-10,222-10;...   %greenish blue. sucrose
                 220-10, 225-10, 248-10;...  %purple. sucrose+ouab  
                 172+5,221+10,222+10;...   %lighter greenish blue. Na glutamate
                 220+5, 225+5, 250] ;  %lighter purple. Na glutamate+ouab               

COLORS = 1/255 * [  95  120  202 ; ...
    93  187 70  ; ...
    241 156 31  ; ...
    237 28  36  ; ...
    129 41  134 ; ...
    75 186 233 ; ...
    102, 141, 60; ... %green
    188,154,125;... %brown/pink
    125,188,186;  ... %complementary to green
    125, 127, 188]; %complementary to brown
 Grey = [0.4 0.4 0.4];
%%
nao=128E-3;
ko=4E-3;
clo=nao+ko;
so=0;
V=-48E-3;
x=2.6E-14; % moles of impermeant intracellular molecules (concentration will depend on volume)
xo=x/50;
z=-1; %average charge of impermeants
wss=cellVolFn(nao,ko,clo,so,V,x,z); %% normal condition. Define normal volume of cell in Liters.
fss=0.7;  
woss=wss/fss-wss;
wtot=wss+woss; %% total volum in liters. This will remail constant;

Do=1; %/mu m^2 / ms intrinsic diffusivity in the ECS
Db=1; %/mu m^2 / ms diffusivity of the intracellular compartment not restricted
Dc=0.1; %/mu m^2 / ms diffusivity of the intracellular compartment restricted
% AXRt=.3;%100/1000; %1/ms
% AXRg=.03;
%Pfactor=0.001;

%Ea_t=52; %47;  %50% kJ/mol
%Ea_g=22;
%A_t=exp(26.62) ; %24 % 26$[1/s]
%A_g=exp(12.52);
%T=25+273.15;
%AXRt=AXR_T(A_t,Ea_t,T)*1E-3; %[1/ms]
%AXRg=AXR_T(A_g,Ea_g,T)*1E-3; %[1/ms]
AXRt=300*1E-3; %[1/ms]
AXRg=30*1E-3;



%% sucrose addition to normal sample
V=-48E-3;
so=0;
    w=cellVolFnPfnf(nao,ko,clo,so,V,wtot,xo,x,z);
    fibo_=w/wtot;
    fa=1-fibo_;
    fb=(1-fa)/2; %% 1/2 and 1/2
    fc=fb;
    Da=Do;
    ADCbo=fa*Da+fb*Db+fc*Dc;
    SxON=multisitekfn(Da,Db,Dc,AXRt,AXRg,fa,fb,fc,1);
    %kbo_tort2(i)=multisitekfn(Da_tort2,Db,Dc,AXRt,AXRg,fa,fb,fc);
    eig1ON=AXRg*(fb+fc)+AXRt*fa;
    eig2ON=AXRt;
    feig1ON=fb+fc;
    feig2ON=fa;
so=.10;
    w=cellVolFnPfnf(nao,ko,clo,so,V,wtot,xo,x,z);
    fibo_=w/wtot;
    fa=1-fibo_;
    fb=(1-fa)/2; %% 1/2 and 1/2
    fc=fb;
    Da=Do;
    ADCbo=fa*Da+fb*Db+fc*Dc;
    SxON100=multisitekfn(Da,Db,Dc,AXRt,AXRg,fa,fb,fc,1);
    %kbo_tort2(i)=multisitekfn(Da_tort2,Db,Dc,AXRt,AXRg,fa,fb,fc);
    eig1ON100=AXRg*(fb+fc)+AXRt*fa;
    eig2ON100=AXRt;
    feig1ON100=fb+fc;
    feig2ON100=fa;
%% sucrose addition to ouabain-treated sample

V=-10E-3;
so=0;
    w=cellVolFnPfnf(nao,ko,clo,so,V,wtot,xo,x,z);
    fiao_=w/wtot;
    fa=1-fiao_;
    fb=(1-fa)/2; %% 1/2 and 1/2
    fc=fb;
    Da=Do;
    ADCao=fa*Da+fb*Db+fc*Dc;
    SxOFF=multisitekfn(Da,Db,Dc,AXRt,AXRg,fa,fb,fc,1);
   % kao_tort2(i)=multisitekfn(Da_tort2,Db,Dc,AXRt,AXRg,fa,fb,fc);
    eig1OFF=AXRg*(fb+fc)+AXRt*fa;
    eig2OFF=AXRt;
    feig1OFF=fb+fc;
    feig2OFF=fa;
%%



%% initialize fig

fig                 = figure();
fig.Units           = 'centimeters';
fig.PaperUnits      = 'centimeters';
fig.Position        = [0 0 8 9.5];
fig.PaperPosition   = fig.Position;

FontName            = 'helvetica';
FontSize            = 7;
FontWeight          = 'normal';



%% PUMP OFF vs PUMP ON vs PUMP ON + 100 mOsm
%% plot I_3 minus T1 fit for individual representative sample
h=subplot(2,1,1);
h.Units             = 'centimeters';
hold on
h.FontName=FontName;
h.FontSize=FontSize;
%h.ylim=ylimhold

h.YLabel.String='DEXR signals ';
h.XLabel.String='$t_m$ [ms]';
%h.XLim=[0 ttotal]; 
h.Position(4)          = h.Position(4)+0.5;
%h.YScale='log';

%%% annotation
ha=annotation('textbox');
ha.Interpreter='latex';
ha.String='A';
ha.FontSize=11;
ha.Units='centimeters';
ha.Position(1)=h.Position(1)-1.1;
ha.Position(2)=h.Position(2)+3.5;
ha.Position(3)=0.5;
ha.Position(4)=0.5;
ha.LineStyle='none'
% ha.Color=[0 0 0];
% ha.BackgroundColor=[.95 .95 .95];
% ha.EdgeColor= [.8 0 0];


%%% choosing representative sample
si=1;
cii=1;
set=1;

y=SxOFF(si).T(cii).I_AXR_I3_subbi(set,:);
x=SxON(si).MixingTimelist;
hl1=plot(x,y);
hl1.LineStyle='none';
hl1.Marker='o';
hl1.Color=COLORS(3,:);
hl1.DisplayName='$V=-10$ mV';

yfit=SxOFF(si).T(cii).fit_I3_subbi{set,1}.Imodel;
hl=plot(x,yfit);
hl.Color=COLORS(3,:)*.75;;



y=SxON(si).T(cii).I_AXR_I3_subbi(set,:);
hl2=plot(x,y);
hl2.LineStyle='none';
hl2.Marker='s';
hl2.Color=COLORS(7,:);
hl2.DisplayName='$V=-48$ mV';

yfit=SxON(si).T(cii).fit_I3_subbi{set,1}.Imodel;
hl=plot(x,yfit);
hl.Color=COLORS(7,:)*.5;

y=SxON100(si).T(cii).I_AXR_I3_subbi(set,:);
hl3=plot(x,y);
hl3.LineStyle='none';
hl3.Marker='^';
hl3.Color=COLORS(1,:);
hl3.DisplayName='$V=-48$ mV + 100 mOsm';

yfit=SxON100(si).T(cii).fit_I3_subbi{set,1}.Imodel;
hl=plot(x,yfit);
hl.Color=COLORS(1,:)*.5;

biexpmodel=0.8*exp(-AXRt*x)+0.2*exp(-AXRg*x);
%biexpmodel=0.7*exp(-eig1ON*x)+0.3*exp(-eig2ON*x);
biexpmodel_=biexpmodel*SxON.T.I0_I3_subbi+SxON.T.baseline_I3_subbi;

% hl=plot(x,biexpmodel_);
% hl.Color=[0 0 0];
% hl.LineStyle='--';


h.YLim(2)=0.33;
legend([hl1,hl2,hl3]);
h.TickDir='Out';

%% plot I_3 minus T1 fit for individual representative sample, log
h=subplot(2,1,2);
h.Units             = 'centimeters';
hold on
h.FontName=FontName;
h.FontSize=FontSize;
%ylimhold=[0 200];

h.YLabel.String='DEXR signals ';
h.XLabel.String='$t_m$ [ms]';
%h.XLim=[0 ttotal]; 
h.Position(4)          = h.Position(4)+0.5;
h.YScale='log';

%%% annotation
ha=annotation('textbox');
ha.Interpreter='latex';
ha.String='B';
ha.FontSize=11;
ha.Units='centimeters';
ha.Position(1)=h.Position(1)-1.1;
ha.Position(2)=h.Position(2)+3.5;
ha.Position(3)=0.5;
ha.Position(4)=0.5;
ha.LineStyle='none'
% ha.Color=[0 0 0];
% ha.BackgroundColor=[.95 .95 .95];
% ha.EdgeColor= [.8 0 0];


%%% choosing representative sample
bl=SxOFF(si).T(cii).fit_I3_subbi{set,1}.baseline; 
I0=SxOFF(si).T(cii).fit_I3_subbi{set,1}.I0

y=(SxOFF(si).T(cii).I_AXR_I3_subbi(set,:)-bl)/I0;
hl=plot(x,y);
hl.LineStyle='none';
hl.Marker='o';
hl.Color=COLORS(3,:);

yfit=(SxOFF(si).T(cii).fit_I3_subbi{set,1}.Imodel-bl)/I0;
hl=plot(x,yfit);
hl.Color=COLORS(3,:)*.75;;

biexpmodel=feig1OFF*exp(-eig1OFF*x)+feig2OFF*exp(-eig2OFF*x);
hl=plot(x,biexpmodel*y(1));
hl.Color=[0 0 0];
hl.LineStyle='--';

bl=SxON(si).T(cii).fit_I3_subbi{set,1}.baseline; 
I0=SxON(si).T(cii).fit_I3_subbi{set,1}.I0

y=(SxON(si).T(cii).I_AXR_I3_subbi(set,:)-bl)/I0;
hl=plot(x,y);
hl.LineStyle='none';
hl.Marker='s';
hl.Color=COLORS(7,:);

yfit=(SxON(si).T(cii).fit_I3_subbi{set,1}.Imodel-bl)/I0;
hl=plot(x,yfit);
hl.Color=COLORS(7,:)*.5;

biexpmodel=feig1ON*exp(-eig1ON*x)+feig2ON*exp(-eig2ON*x);
hl=plot(x,biexpmodel*y(1));
hl.Color=[0 0 0];
hl.LineStyle='--';

bl=SxON100(si).T(cii).fit_I3_subbi{set,1}.baseline; 
I0=SxON100(si).T(cii).fit_I3_subbi{set,1}.I0

y=(SxON100(si).T(cii).I_AXR_I3_subbi(set,:)-bl)/I0;
hl=plot(x,y);
hl.LineStyle='none';
hl.Marker='^';
hl.Color=COLORS(1,:);

yfit=(SxON100(si).T(cii).fit_I3_subbi{set,1}.Imodel-bl)/I0;
hl=plot(x,yfit);
hl.Color=COLORS(1,:)*.5;

biexpmodel=feig1ON100*exp(-eig1ON100*x)+feig2ON100*exp(-eig2ON100*x);
hl=plot(x,biexpmodel*y(1));
hl.Color=[0 0 0];
hl.LineStyle='--';

h.TickDir='Out'
h.YLim(1)=0.01;
h.XLim(2)=50; %50
return
fig.Renderer='Painters';
print(fig,'figure11.tif','-dtiff','-r300')
print(fig,'figure11.eps','-depsc')
