%% Fig 3 source code
%% Williamson, et al, Magnetic Resonance Letters (2025)
%%%  figure showing represenattive raw exchange data from individual simulations.
%%% requires the folder "functions" to be added to the path
%%% note that outS=1 in multisitekfn() so that it will output the fits

clear all
%close all

set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
%% simulate data
clear all
close all

set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
    

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
%% Define ion concentration, moles and charge of impermeant osmolytes, and voltage
%%% these are defined and will remail constant
nao=128E-3; %intracellular 
ko=4E-3;
clo=nao+ko;
so=0;
V=-48E-3;  % voltage 
x=2.6E-14; % moles of impermeant intracellular molecules (concentration will depend on volume)
xo=x/50;   % moles of impermeant extracellular molecules
z=-1; %average charge of impermeants

%%% wss and fss are initial conditions needed to set wtot. They will change. 
wss=cellVolFn(nao,ko,clo,so,V,x,z); % find initial normal volume of cell in Liters. This uses the pump
fss=0.7;   % Define initial normal intracellular fraction.
woss=wss/fss-wss; % find initial normal intracellular fraction.
%%%
wtot=wss+woss; %% total volume in liters. This is defined and will remail constant;

Da=1; %/mu m^2 / ms  diffusivity in the ECS
Db=1; %/mu m^2 / ms diffusivity of the intracellular compartment not restricted
Dc=0.1; %/mu m^2 / ms diffusivity of the intracellular compartment restricted

AXRt=300*1E-3; %[1/ms]  %%% used to define transmembrane exchange kinetics in exchange matrix K
AXRg=30*1E-3; %[1/ms]  %%% used to define transmembrane exchange kinetics in exchange matrix K


so_=0;
%% sucrose addition to normal sample
V=-48E-3;
    so=so_
    w=cellVolFnPfnf(nao,ko,clo,so,V,wtot,xo,x,z);
    fibo_=w/wtot;
    fa=1-fibo_;
    fb=(1-fa)/2; %% 1/2 and 1/2
    fc=fb;
    ADCbo=fa*Da+fb*Db+fc*Dc;
    SxON=multisitekfn(Da,Db,Dc,AXRt,AXRg,fa,fb,fc,1);


%% initialize fig

fig                 = figure();
fig.Units           = 'centimeters';
fig.PaperUnits      = 'centimeters';
fig.Position        = [0 0 8 8];
fig.PaperPosition   = fig.Position;

FontName            = 'helvetica';
FontSize            = 7;
FontWeight          = 'normal';



%% PUMP OFF vs PUMP ON vs PUMP ON + 100 mOsm
%% plot I_2 for an individual representative sample
h=subplot(1,1,1);
h.Units             = 'centimeters';
hold on
h.FontName=FontName;
h.FontSize=FontSize;
ylimhold=[0.05 0.261];
h.YLim=ylimhold;
%h.YScale='log';

h.YLabel.String='simulated signals';
h.XLabel.String='$t_m$ [ms]';
%h.XLim=[0 ttotal]; 
h%.Position(4)          = h.Position(4)+0.5;
%h.Title.String= 'PUMP ON vs. PUMP ON + 100 mOsm';
h.Title.FontWeight='normal';

%%% choosing representative sample
si=1;  %% Sample index fixed
cii=1; %% 25C
set=1;


y=SxON(si).T(cii).I_2(set,:);
x=SxON(si).MixingTimelist;

hl=plot(x,y);
hl.LineStyle='none';
hl.Marker='<';
%hl.MarkerFaceColor=COLORS(2,:);
hl.Color=COLORS(6,:);
hl.DisplayName='$I(0.200, 0.735)$ ';

yfit=SxON(si).T(cii).fit_I2biexp{set,1}.Imodel;
hl=plot(x,yfit);
hl.Color=COLORS(6,:);
hl.DisplayName='$I(t_m)=I_0( w_1\exp{(-t_m R_1^1)} + (1-w_1)\exp{(-t_m R_1^2)}$'

y=SxON(si).T(cii).I_3(set,:);
hl=plot(x,y);
hl.LineStyle='none';
hl.Marker='>';
hl.Color=COLORS(5,:);
hl.DisplayName='$I(0.593, 0.580)$';

y=SxON(si).T(cii).I_AXR_I3_subbi(set,:);
hl=plot(x,y);
hl.LineStyle='none';
hl.Marker='s';
hl.Color=COLORS(7,:);
hl.DisplayName='DEXR signal';

yfit=SxON(si).T(cii).fit_I3_subbi{set,1}.Imodel;
hl=plot(x,yfit);
hl.Color=COLORS(7,:)*.5;
hl.DisplayName='exchange model';
h.YLim=[0 .25]
hl.DisplayName='$I(t_m)=I_0\exp{(-t_m k)}+B$'

legend('Location','South')

%fig.Renderer='Painters';
return
%% move legend before saving
print(fig,'figure3.tif','-dtiff','-r300')
print(fig,'figure3.eps','-depsc')
