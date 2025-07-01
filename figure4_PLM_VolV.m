%% Fig 4 source code
%% Williamson, et al, Magnetic Resonance Letters (2025)
%%% use Pump Leak model to predict a bar chart of volume and voltage
%%% Pump leak model branches from Kay, Froniters (2017) https://doi.org/10.3389/fcell.2017.00041

clear all
%close all

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

%https://htmlcolorcodes.com/color-picker/             
%COLORSpatchnew=1/255*[163,213,105];


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

%% make arrays with data

%%% define aCSF
nao=128;
ko=4;
clo=128+4;
xo=0; 
[Warray(1),Varray(1),Warray(2),Varray(2)]=PLEfn(nao*1E-3,ko*1E-3,clo*1E-3,xo*1E-3);
%%% define aCSF+100 mOsm
xo=100; 
[Warray(3),Varray(3),Warray(4),Varray(4)]=PLEfn(nao*1E-3,ko*1E-3,clo*1E-3,xo*1E-3);
%%% define 0 NaCl +257 mM sucrose
nao=.1;
ko=4;
clo=.1+4;
xo=256; 
[Warray(5),Varray(5),Warray(6),Varray(6)]=PLEfn(nao*1E-3,ko*1E-3,clo*1E-3,xo*1E-3);
%%% define 0 NaCl + 128.35 mM Na gluconate 
nao=128;
ko=4;
clo=0+4;
xo=128; 
[Warray(7),Varray(7),Warray(8),Varray(8)]=PLEfn(nao*1E-3,ko*1E-3,clo*1E-3,xo*1E-3);

%%% normalizing volume
Warray=Warray/Warray(1)*100
%%

set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');


fig                 = figure();
fig.Units           = 'centimeters';
fig.PaperUnits      = 'centimeters';
fig.Position        = [0 0 18 11];

fig.PaperPosition   = fig.Position;

FontName            = 'helvetica';
FontSize            = 7;
FontWeight          = 'normal';

%% cell volume bar chart
h=subplot(1,2,1);
h.Units='centimeters';
h.FontName=FontName;
h.FontSize=FontSize;
h.FontName=FontName;
hold on

h.YLabel.String='cell volume / normal cell volume x 100 [\%]';
h.XLabel.String='';
h.XTick=[1:8];
h.XTickLabel={'normal 128 mM NaCl, PUMP ON','128 mM NaCl, PUMP OFF', '128 mM NaCl + 100 mM osmolyte, PUMP ON', '128 mM NaCl + 100 mM osmolyte, PUMP OFF ',' 0 NaCl, 256 mM osmolyte, PUMP ON',' 0 NaCl, 256 mM osmolyte, PUMP OFF',' 0 NaCl, 128 mM sodium gluconate, PUMP ON', '0 NaCl, 128 mM sodium gluconate, PUMP OFF'};

%h.XLim=[0.5 4.5];
h.YLim=[0 200];
%  h.Position(3)=h.Position(3)+1.2
  h.Position(1)=h.Position(1)+.2
 h.Position(2)=h.Position(2)+2;
 h.Position(4)=h.Position(4)-2
% h.Position(4)=h.Position(4)-0.1;
ha=annotation('textbox');
ha.Units='centimeters';
ha.Interpreter='latex';
    ha.String='A';
    ha.Position(1)=h.Position(1)-1.2;
    ha.Position(2)=h.Position(2)+7;
ha.FontSize=11;
%ha.Position(2)=h.Position(2)+6.1;%4.15;
ha.Position(3)=0.5;
ha.Position(4)=0.5;
ha.Color=[0 0 0];
%ha.BackgroundColor=[.95 .95 .95];
ha.EdgeColor= 'none';% [.8 0 0];


xtickangle(30);

x=[1:8];
y=Warray;

hb=bar(x,y);
hb.FaceColor = 'flat';
%hb.FaceAlpha = .4;
 hb.CData(1,:) = COLORSpatch(1,:);
hb.CData(2,:) = COLORSpatch(3,:);
hb.CData(3,:) = COLORSpatch(5,:);
hb.CData(4,:) = COLORSpatch(6,:);
hb.CData(5,:) = COLORSpatch(7,:);
hb.CData(6,:) = COLORSpatch(8,:);
hb.CData(7,:) = COLORSpatch(9,:);
hb.CData(8,:) = COLORSpatch(10,:);
%% voltage bar chart
h=subplot(1,2,2);
h.Units='centimeters';
h.FontName=FontName;
h.FontSize=FontSize;
h.FontName=FontName;
hold on
h.YLabel.String='voltage [mV]';
h.XLabel.String='';
h.XTick=[1:8];
h.XTickLabel={'normal 128 mM NaCl, PUMP ON','128 mM NaCl, PUMP OFF', '128 mM NaCl + 100 mM osmolyte, PUMP ON', '128 mM NaCl + 100 mM osmolyte, PUMP OFF ',' 0 NaCl, 256 mM osmolyte, PUMP ON',' 0 NaCl, 256 mM osmolyte, PUMP OFF',' 0 NaCl, 128 mM sodium gluconate, PUMP ON', '0 NaCl, 128 mM sodium gluconate, PUMP OFF'};

%h.XLim=[0.5 4.5];
%.YLim=[0 180];
%  h.Position(3)=h.Position(3)+1.2
%  h.Position(1)=h.Position(1)+.2
 h.Position(2)=h.Position(2)+2;
 h.Position(4)=h.Position(4)-2
% h.Position(4)=h.Position(4)-0.1;
ha=annotation('textbox');
ha.Units='centimeters';
ha.Interpreter='latex';
    ha.String='B';
    ha.Position(1)=h.Position(1)-1.2;
    ha.Position(2)=h.Position(2)+7;
ha.FontSize=11;
%ha.Position(2)=h.Position(2)+6.1;%4.15;
ha.Position(3)=0.5;
ha.Position(4)=0.5;
ha.Color=[0 0 0];
%ha.BackgroundColor=[.95 .95 .95];
ha.EdgeColor= 'none';% [.8 0 0];


xtickangle(30);

x=[1:8];
y=Varray;

hb=bar(x,y);
hb.FaceColor = 'flat';
%hb.FaceAlpha = .4;
 hb.CData(1,:) = COLORSpatch(1,:);
hb.CData(2,:) = COLORSpatch(3,:);
hb.CData(3,:) = COLORSpatch(5,:);
hb.CData(4,:) = COLORSpatch(6,:);
hb.CData(5,:) = COLORSpatch(7,:);
hb.CData(6,:) = COLORSpatch(8,:);
hb.CData(7,:) = COLORSpatch(9,:);
hb.CData(8,:) = COLORSpatch(10,:);

%%
print(fig,'figure4.eps','-depsc')
print(fig,'figure4.png','-dpng','-r300')
