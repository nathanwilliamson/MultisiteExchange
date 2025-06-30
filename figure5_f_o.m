%% Fig 5 source code
%% Williamson, et al, Magnetic Resonance Letters (2025)

%%% Predict fo and osmolarities of trapped impermeants for a steady state model of cell volume with an added pressure from extracellular osmolytes. 
%%% assumes equal osmolarity and electroneutrality  


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

AXRt=300*1E-3; %[1/ms]
AXRg=30*1E-3;


so_=linspace(0,300E-3,51);
so_=horzcat(so_,linspace(301E-3,2000E-3,49))
%% sucrose addition to normal sample
V=-48E-3;
for i=1:length(so_)
    so=so_(i);
    w=cellVolFnPfnf(nao,ko,clo,so,V,wtot,xo,x,z);
    fibo_(i)=w/wtot;
    fa=1-fibo_(i);
    fb=(1-fa)/2; %% 1/2 and 1/2
    fc=fb;
    Da=Do;

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
    Da=Do;

end

%% make plots

set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');


fig                 = figure();
fig.Units           = 'centimeters';
fig.PaperUnits      = 'centimeters';
fig.Position        = [0 0 18 5.9];

fig.PaperPosition   = fig.Position;

FontName            = 'helvetica';
FontSize            = 7;
FontWeight          = 'normal';

%% fo
h=subplot(1,3,1);
h.Units='centimeters';
h.FontName=FontName;
h.FontSize=FontSize;
h.FontName=FontName;
hold on

h.XLabel.String='+ osmolyte [mOsm]';
h.YLabel.String='$f_\mathrm{o}$';

hl=plot(so_(1:51)*1000,1-fiao_(1:51))
hl.Color=COLORSpatch(6,:) -[0 50 0]/255;
hl.LineWidth=1.5;

hl=plot(so_(1:51)*1000,1-fibo_(1:51))
hl.Color=COLORS(7,:);


ha=annotation('textbox');
ha.Units='centimeters';
ha.Interpreter='latex';
    ha.String='A';
    ha.Position(1)=h.Position(1)-1;
    ha.Position(2)=h.Position(2)+4.5;
ha.FontSize=11;
%ha.Position(2)=h.Position(2)+6.1;%4.15;
ha.Position(3)=0.5;
ha.Position(4)=0.5;
ha.Color=[0 0 0];
%ha.BackgroundColor=[.95 .95 .95];
ha.EdgeColor= 'none';% [.8 0 0];


%% osmolarity of extracellular impermeants

h=subplot(1,3,2);
h.Units='centimeters';
h.FontName=FontName;
h.FontSize=FontSize;
h.FontName=FontName;
hold on


h.XLabel.String='+osmolyte [mOsm]';
h.YLabel.String     = '$x_\mathrm{o}/(f_\mathrm{o}w_\mathrm{tot})$ [mOsm]';

hl=plot(so_(1:49)*1000,xo./((1-fiao_(1:49))*wtot)*1000)
hl.Color=COLORSpatch(6,:) -[0 50 0]/255;
hl.LineWidth=1.5;

hl=plot(so_(1:49)*1000,xo./((1-fibo_(1:49))*wtot)*1000)
hl.Color=COLORS(7,:)
% hl=plot(So_*1000,kbo_tort2)
% hl.Color=COLORS(7,:)/1.6;
% hl.LineStyle='--';  

 %P = polyfit(ADCao,kao_tort1*1000,1);
 %f = polyval(P,ADCao);
 % plot(ADCao,f,'--')


ha=annotation('textbox');
ha.Units='centimeters';
ha.Interpreter='latex';
    ha.String='B';
    ha.Position(1)=h.Position(1)-1;
    ha.Position(2)=h.Position(2)+4.5;
ha.FontSize=11;
%ha.Position(2)=h.Position(2)+6.1;%4.15;
ha.Position(3)=0.5;
ha.Position(4)=0.5;
ha.Color=[0 0 0];
%ha.BackgroundColor=[.95 .95 .95];
ha.EdgeColor= 'none';% [.8 0 0];
%% osmolarity of intracellular impermeants

h=subplot(1,3,3);
h.Units='centimeters';
h.FontName=FontName;
h.FontSize=FontSize;
h.FontName=FontName;
hold on


h.XLabel.String='+osmolyte [mOsm]';
h.YLabel.String     = '$x_\mathrm{i}/(f_\mathrm{i}w_\mathrm{tot})$ [mOsm]';

hl=plot(so_(1:49)*1000,x./(fibo_(1:49)*wtot)*1000);
hl.Color=COLORS(7,:)

hl=plot(so_(1:49)*1000,x./(fiao_(1:49)*wtot)*1000);
hl.Color=COLORSpatch(6,:) -[0 50 0]/255;
hl.LineWidth=1.5;


% hl=plot(So_*1000,kbo_tort2)
% hl.Color=COLORS(7,:)/1.6;
% hl.LineStyle='--';  
 legend('$V=-48$ mV', '$V=-10$ mV',Location='northwest')

 %P = polyfit(ADCao,kao_tort1*1000,1);
 %f = polyval(P,ADCao);
 % plot(ADCao,f,'--')

ha=annotation('textbox');
ha.Units='centimeters';
ha.Interpreter='latex';
    ha.String='C';
    ha.Position(1)=h.Position(1)-1;
    ha.Position(2)=h.Position(2)+4.5;
ha.FontSize=11;
%ha.Position(2)=h.Position(2)+6.1;%4.15;
ha.Position(3)=0.5;
ha.Position(4)=0.5;
ha.Color=[0 0 0];
%ha.BackgroundColor=[.95 .95 .95];
ha.EdgeColor= 'none';% [.8 0 0];

%% total c
%% osmolarity of intracellular impermeants
% 
% h=subplot(3,3,3);
% h.Units='centimeters';
% h.FontName=FontName;
% h.FontSize=FontSize;
% h.FontName=FontName;
% hold on
% 
% 
% h.XLabel.String='+Osmolyte [mOsm]';
% h.YLabel.String     = 'c [mM]';
% 
% hl=plot(so_(1:49)*1000,(xo./((1-fiao_(1:49))*wtot)+so_(1:49)+nao+ko+clo)*1000)
% hl.Color=COLORSpatch(6,:) -[0 50 0]/255;
% hl.LineWidth=1.5;
% 
% hl=plot(so_(1:49)*1000,(xo./((1-fibo_(1:49))*wtot)+so_(1:49)+nao+ko+clo)*1000)
% hl.Color=COLORS(7,:)
% 
% % hl=plot(So_*1000,kbo_tort2)
% % hl.Color=COLORS(7,:)/1.6;
% % hl.LineStyle='--';  
% 
% 
%  %P = polyfit(ADCao,kao_tort1*1000,1);
%  %f = polyval(P,ADCao);
%  % plot(ADCao,f,'--')
% 
% ha=annotation('textbox');
% ha.Units='centimeters';
% ha.Interpreter='latex';
%     ha.String='C';
%     ha.Position(1)=h.Position(1)-1;
%     ha.Position(2)=h.Position(2)+4.3;
% ha.FontSize=11;
% %ha.Position(2)=h.Position(2)+6.1;%4.15;
% ha.Position(3)=0.5;
% ha.Position(4)=0.5;
% ha.Color=[0 0 0];
% %ha.BackgroundColor=[.95 .95 .95];
% ha.EdgeColor= 'none';% [.8 0 0];



%%


print(fig,'figure5.eps','-depsc')
print(fig,'figure5.png','-dpng','-r300')

%%
figure 
hold on
h=subplot(1,1,1);
h.Units='centimeters';
h.FontName=FontName;
h.FontSize=FontSize;
h.FontName=FontName;
hold on


h.XLabel.String='+osmolyte [mOsm]';
h.YLabel.String     = 'so + xo/(fo*wtot) + (Nao+Ko+Clo) - (Nai+Ki+Cli) [mM]';

R=8.314; F=96485; % R (RT/F) in Volts, Faraday's constant C/mol
T=273+25;

V=-0.01;
hl=plot(so_(1:49)*1000,(xo./((1-fiao_(1:49))*wtot)+so_(1:49)+2*clo-2*clo*exp(F*V/(R*T)))*1000)
hl.Color=COLORSpatch(6,:) -[0 50 0]/255;
hl.LineWidth=1.5;

V=-0.048
hl=plot(so_(1:49)*1000,(xo./((1-fibo_(1:49))*wtot)+so_(1:49)+nao+ko+clo-2*clo*exp(F*V/(R*T)))*1000)
hl.Color=COLORS(7,:)

