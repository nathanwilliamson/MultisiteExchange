%% Fig 6 source code
%% Williamson, et al, Magnetic Resonance Letters (2025)
%%% predicts the dependence of ADC and AXR on f_o and osmolarity 
%%% for the case of two-site exchange between sites A and C
clear all
close all

newrun_yn=1; %% set to 1 if running for the first time. 
            %%%set to 0 if already ran and variables are saved in workspace

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

Da=1; %/mu m^2 / ms intrinsic diffusivity in the ECS
Db=1; %/mu m^2 / ms diffusivity of the intracellular compartment not restricted
Dc=0.1; %/mu m^2 / ms diffusivity of the intracellular compartment restricted

AXRt=300*1E-3; %[1/ms]
AXRg=30*1E-3;

if newrun_yn==1
so_=linspace(0,300E-3,51);
so_=horzcat(so_,linspace(301E-3,2000E-3,49))
%% sucrose addition to normal sample
V=-48E-3;
for i=1:length(so_)
    so=so_(i);
    w=cellVolFnPfnf(nao,ko,clo,so,V,wtot,xo,x,z);
    fibo_(i)=w/wtot;
    fa=1-fibo_(i);
    fc=(1-fa); 
    ADCbo(i)=fa*Da+fc*Dc;
    kbo(i)=twositekfn(Da,Dc,AXRt,fc);
end
%% sucrose addition to ouabain-treated sample

V=-10E-3;
for i=1:length(so_)
    so=so_(i);
    w=cellVolFnPfnf(nao,ko,clo,so,V,wtot,xo,x,z);
    fiao_(i)=w/wtot;
    fa=1-fiao_(i);
    fc=(1-fa); %% 1/2 and 1/2
    ADCao(i)=fa*Da+fc*Dc;
    kao(i)=twositekfn(Da,Dc,AXRt,fc);
end
%% save data
save('2siteModel_baseCase')
elseif newrun_yn==0
    load('2siteModel_baseCase.mat')
end
%% make plots

set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');


fig                 = figure();
fig.Units           = 'centimeters';
fig.PaperUnits      = 'centimeters';
fig.Position        = [0 0 18 13];

fig.PaperPosition   = fig.Position;

FontName            = 'helvetica';
FontSize            = 7;
FontWeight          = 'normal';


%% ADC vs fo
h=subplot(2,2,1);
h.Units='centimeters';
h.FontName=FontName;
h.FontSize=FontSize;
h.FontName=FontName;
hold on

h.XLabel.String='$f_\mathrm{o}$';
h.YLabel.String='ADC [$\mathrm{\mu m^2/ms}$]';

% hl=plot(1-fibo_(1:49),ADCbo(1:49));
% hl.Color=COLORS(7,:);

hl=plot(1-fiao_(1:49),ADCao(1:49))
hl.Color='k';
%hl.LineWidth=1.5;
hl=plot(1-fiao_(52:end),ADCao(52:end))
hl.Color=[.5 .5 .5];
hl.LineStyle='--';

ha=annotation('textbox');
ha.Units='centimeters';
ha.Interpreter='latex';
    ha.String='A';
    ha.Position(1)=h.Position(1)-1;
    ha.Position(2)=h.Position(2)+4.3;
ha.FontSize=11;
%ha.Position(2)=h.Position(2)+6.1;%4.15;
ha.Position(3)=0.5;
ha.Position(4)=0.5;
ha.Color=[0 0 0];
%ha.BackgroundColor=[.95 .95 .95];
ha.EdgeColor= 'none';% [.8 0 0];

legend('0--300 mOsm','300--2000 mOsm','Location','southeast')
%% k vs fo
h=subplot(2,2,2);
h.Units='centimeters';
h.FontName=FontName;
h.FontSize=FontSize;
h.FontName=FontName;
hold on

h.XLabel.String='$f_\mathrm{o}$';
h.YLabel.String= 'AXR [$\mathrm{s}^{-1}$]';

hl=plot(1-fiao_(1:49),kao(1:49)*1000)
hl.Color='k';
%hl.LineWidth=1.5;
hl=plot(1-fiao_(52:end),kao(52:end)*1000)
hl.Color=[.5 .5 .5];
hl.LineStyle='--';


kgroundtruth=AXRt.*fiao_+AXRt.*(1-fiao_); %ground truth for two site exchange model

hl=plot((1-fiao_), kgroundtruth*1000);
hl.Color='r';
hl.LineStyle=':';

h.YLim=[0 320];

legend('0--300 mOsm','300--2000 mOsm','ground truth $k$','Location','southeast')

ha=annotation('textbox');
ha.Units='centimeters';
ha.Interpreter='latex';
    ha.String='B';
    ha.Position(1)=h.Position(1)-1;
    ha.Position(2)=h.Position(2)+4.3;
ha.FontSize=11;
%ha.Position(2)=h.Position(2)+6.1;%4.15;
ha.Position(3)=0.5;
ha.Position(4)=0.5;
ha.Color=[0 0 0];
%ha.BackgroundColor=[.95 .95 .95];
ha.EdgeColor= 'none';% [.8 0 0];
%%
h=subplot(2,3,4);
h.Units='centimeters';
h.FontName=FontName;
h.FontSize=FontSize;
h.FontName=FontName;
hold on

h.XLabel.String='+ osmolyte [mOsm]';
h.YLabel.String='ADC [$\mathrm{\mu m^2/ms}$]';
hl=plot(so_(1:51)*1000,ADCbo(1:51))
hl.Color=COLORS(7,:)

% hl=plot(So_*1000,ADCbo_tort2)
% hl.Color=COLORS(7,:)/1.6;
% hl.LineStyle='--';  

hl=plot(so_(1:51)*1000,ADCao(1:51))
hl.Color=COLORSpatch(6,:) -[0 50 0]/255;
hl.LineWidth=1.5;


ha=annotation('textbox');
ha.Units='centimeters';
ha.Interpreter='latex';
    ha.String='C';
    ha.Position(1)=h.Position(1)-1;
    ha.Position(2)=h.Position(2)+4.3;
ha.FontSize=11;
%ha.Position(2)=h.Position(2)+6.1;%4.15;
ha.Position(3)=0.5;
ha.Position(4)=0.5;
ha.Color=[0 0 0];
%ha.BackgroundColor=[.95 .95 .95];
ha.EdgeColor= 'none';% [.8 0 0];

%%
h=subplot(2,3,5);
h.Units='centimeters';
h.FontName=FontName;
h.FontSize=FontSize;
h.FontName=FontName;
hold on

h.XLabel.String='+ osmolyte [mOsm]';
h.YLabel.String     = 'AXR [$\mathrm{s^{-1}}$]';

hl=plot(so_(1:51)*1000,kao(1:51)*1000)
hl.Color=COLORSpatch(6,:) -[0 50 0]/255;
hl.LineWidth=1.5;

hl=plot(so_(1:51)*1000,kbo(1:51)*1000)
hl.Color=COLORS(7,:)

% hl=plot(So_*1000,kbo_tort2)
% hl.Color=COLORS(7,:)/1.6;
% hl.LineStyle='--';  

ha=annotation('textbox');
ha.Units='centimeters';
ha.Interpreter='latex';
    ha.String='D';
    ha.Position(1)=h.Position(1)-1;
    ha.Position(2)=h.Position(2)+4.3;
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

h.YLim=[0 320];
%% correlation ADC and k

h=subplot(2,3,6);
h.Units='centimeters';
h.FontName=FontName;
h.FontSize=FontSize;
h.FontName=FontName;
hold on


h.XLabel.String='ADC [$\mathrm{\mu m^2/ms}$]';
h.YLabel.String     = 'AXR [$\mathrm{s}^{-1}$]';
%%%
hl2=plot(ADCao(1:51),kao(1:51)*1000)
hl2.Color=COLORSpatch(6,:) -[0 50 0]/255;
hl2.LineWidth=1.5;

hl1=plot(ADCbo(1:51),kbo(1:51)*1000)
hl1.Color=COLORS(7,:)


% hl=plot(So_*1000,kbo_tort2)
% hl.Color=COLORS(7,:)/1.6;
% hl.LineStyle='--';  


 P = polyfit(ADCao(1:51),kao(1:51)*1000,1);
 f = polyval(P,ADCao(1:51));
 % plot(ADCao,f,'--')

ha=annotation('textbox');
ha.Units='centimeters';
ha.Interpreter='latex';
    ha.String='E';
    ha.Position(1)=h.Position(1)-1;
    ha.Position(2)=h.Position(2)+4.3;
ha.FontSize=11;
%ha.Position(2)=h.Position(2)+6.1;%4.15;
ha.Position(3)=0.5;
ha.Position(4)=0.5;
ha.Color=[0 0 0];
%ha.BackgroundColor=[.95 .95 .95];
ha.EdgeColor= 'none';% [.8 0 0];
h.YLim=[0 320];
legend('$V=-48$ mV','$V=-10$ mV',Location='south')
legend([hl1,hl2], '$V=-48$ mV','$V=-10$ mV',Location='northwest')

%%

print(fig,'figure6.eps','-depsc')
print(fig,'figure6.png','-dpng','-r300')
