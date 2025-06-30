%% Fig 8, 9, and 10 source code
%% Williamson, et al, Magnetic Resonance Letters (2025)
%%% predicts the dependence of ADC and AXR on f_o and osmolarity 
%%% for additional cases of three-site exchange between sites A, B and C


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


%% need to change Do, Da and Db in the code below 
%%%%% by commenting out/uncommenting code
Do=1.7; % figure 10 %/mu m^2 / ms intrinsic diffusivity  in the ECS (see also lines
%Db=1.5; % Figure 8 %/mu m^2 / ms diffusivity of the intracellular compartment not restricted
%Db=0.5;% Figure 9 %/mu m^2 / ms diffusivity of the intracellular compartment not restricted
 Db=1; %/mu m^2 / ms diffusivity of the intracellular compartment not restricted
%Da=1; %/mu m^2 / ms diffusivity of the extracellular compartment
Dc=0.1; %/mu m^2 / ms diffusivity of the intracellular compartment restricted

AXRt=300*1E-3; %[1/ms]
AXRg=30*1E-3;


so_=linspace(0,300E-3,51);
%% sucrose addition to normal sample
V=-48E-3;
for i=1:length(so_)
    so=so_(i);
    w=cellVolFnPfnf(nao,ko,clo,so,V,wtot,xo,x,z);
    fibo_(i)=w/wtot;
    fa=1-fibo_(i);
    fb=(1-fa)/2; %% 1/2 and 1/2
    fc=fb;
    Da=Do*fa; %Figure 10
    ADCbo(i)=fa*Da+fb*Db+fc*Dc;
    kbo(i)=multisitekfn(Da,Db,Dc,AXRt,AXRg,fa,fb,fc,0);
end
%% sucrose addition to ouabain-treated sample

V=-10E-3;
for i=1:length(so_)
    so=so_(i);
    w=cellVolFnPfnf(nao,ko,clo,so,V,wtot,xo,x,z);
    fiao_(i)=w/wtot;
    fa=1-fiao_(i);
    fb=(1-fa)/2; %% 1/2 and 1/2
    fc=fb;
    Da=Do*fa; %Figure 10
    ADCao(i)=fa*Da+fb*Db+fc*Dc;
    kao(i)=multisitekfn(Da,Db,Dc,AXRt,AXRg,fa,fb,fc,0);
end

%% make plots

set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');


fig                 = figure();
fig.Units           = 'centimeters';
fig.PaperUnits      = 'centimeters';
fig.Position        = [0 0 18 5.5];

fig.PaperPosition   = fig.Position;

FontName            = 'helvetica';
FontSize            = 7;
FontWeight          = 'normal';


h=subplot(1,3,1);
h.Units='centimeters';
h.FontName=FontName;
h.FontSize=FontSize;
h.FontName=FontName;
h.Position=h.Position+[0 .5 0 -.5]
hold on

h.XLabel.String='+ osmolyte [mOsm]';
h.YLabel.String='ADC [$\mathrm{\mu m^2/ms}$]';
hl=plot(so_*1000,ADCbo)
hl.Color=COLORS(7,:)

% hl=plot(So_*1000,ADCbo_tort2)
% hl.Color=COLORS(7,:)/1.6;
% hl.LineStyle='--';  

hl=plot(so_*1000,ADCao)
hl.Color=COLORSpatch(6,:) -[0 50 0]/255;
hl.LineWidth=1.5;


ha=annotation('textbox');
ha.Units='centimeters';
ha.Interpreter='latex';
    ha.String='A';
    ha.Position(1)=h.Position(1)-1;
    ha.Position(2)=h.Position(2)+4;
ha.FontSize=11;
%ha.Position(2)=h.Position(2)+6.1;%4.15;
ha.Position(3)=0.5;
ha.Position(4)=0.5;
ha.Color=[0 0 0];
%ha.BackgroundColor=[.95 .95 .95];
ha.EdgeColor= 'none';% [.8 0 0];

% hl=plot(So_*1000,ADCao_tort2)
% hl.Color=COLORSpatch(3,:)/3;
% hl.LineStyle='--';  

h=subplot(1,3,2);
h.Units='centimeters';
h.FontName=FontName;
h.FontSize=FontSize;
h.FontName=FontName;
h.Position=h.Position+[0 .5 0 -.5]
hold on

h.XLabel.String='+ osmolyte [mOsm]';
h.YLabel.String     = 'AXR [$\mathrm{s}^{-1}$]';

hl=plot(so_*1000,kbo*1000)
hl.Color=COLORS(7,:)

% hl=plot(So_*1000,kbo_tort2)
% hl.Color=COLORS(7,:)/1.6;
% hl.LineStyle='--';  

hl=plot(so_*1000,kao*1000)
hl.Color=COLORSpatch(6,:) -[0 50 0]/255;
hl.LineWidth=1.5;

ha=annotation('textbox');
ha.Units='centimeters';
ha.Interpreter='latex';
    ha.String='B';
    ha.Position(1)=h.Position(1)-1;
    ha.Position(2)=h.Position(2)+4;
ha.FontSize=11;
%ha.Position(2)=h.Position(2)+6.1;%4.15;
ha.Position(3)=0.5;
ha.Position(4)=0.5;
ha.Color=[0 0 0];
%ha.BackgroundColor=[.95 .95 .95];
ha.EdgeColor= 'none';% [.8 0 0];
% hl=plot(So_*1000,kao_tort2)
% hl.Color=COLORSpatch(3,:)/3;
% hl.LineStyle='--';  

%% correlation ADC and k

h=subplot(1,3,3);
h.Units='centimeters';
h.FontName=FontName;
h.FontSize=FontSize;
h.FontName=FontName;
h.Position=h.Position+[0 .5 0 -.5]
hold on


h.XLabel.String='ADC [$\mathrm{\mu m^2/ms}$]';
h.YLabel.String     = 'AXR [$\mathrm{s}^{-1}$]';

hl2=plot(ADCao,kao*1000)
hl2.Color=COLORSpatch(6,:) -[0 50 0]/255;
hl2.LineWidth=1.5;

hl1=plot(ADCbo,kbo*1000)
hl1.Color=COLORS(7,:)
legend([hl1,hl2], '$V=-48$ mV','$V=-10$ mV',Location='southeast')
% hl=plot(So_*1000,kbo_tort2)
% hl.Color=COLORS(7,:)/1.6;
% hl.LineStyle='--';  


 P = polyfit(ADCao,kao*1000,1);
 f = polyval(P,ADCao);
 % plot(ADCao,f,'--')

ha=annotation('textbox');
ha.Units='centimeters';
ha.Interpreter='latex';
    ha.String='C';
    ha.Position(1)=h.Position(1)-1;
    ha.Position(2)=h.Position(2)+4;
ha.FontSize=11;
%ha.Position(2)=h.Position(2)+6.1;%4.15;
ha.Position(3)=0.5;
ha.Position(4)=0.5;
ha.Color=[0 0 0];
%ha.BackgroundColor=[.95 .95 .95];
ha.EdgeColor= 'none';% [.8 0 0];


print(fig,'figure8_.eps','-depsc')
print(fig,'figure8_.png','-dpng','-r300')


