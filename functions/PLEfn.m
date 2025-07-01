%% This code comes from a version originally publised along with Kay, Frontiers (2017) 
%% PLE_Frontiers - solve the PLEs using charge difference method with simple Euler integration
% the simulation starts with the Na-pump off, which is then turned on 
% at t=ton, the pump is then turned off at t=toff 
% p=pump rate, X = moles of impermeant ion, z = average charge of X 
% all concentrations in M, V in volts, dimensions dm , T=25oC 
% start at area & volume defined by radius (rad, in um ), the volume is allowed to change
% but the area is fixed,  w = calculated cell volume 
% intracellular concentrations = k,na,cl,x
% conductances = gna, gk , gcl 
% extracellular concentrations = ko,nao,clo,
% Oso=extracellular osmolarity, Osi=intracellular osmolarity 
% so= concentration of impermeant extracellular molecule with average charge = zso
% water equilibration is assumed to be instantaneous 
% dt=time step 

% Alan R Kay, Dept Biology, University of Iowa -  alan-kay@uiowa.edu 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% additions made for Williamson, et al, Magnetic Resonance Letters (2025)
%%% Turned into a function to predict volume and voltage at steady-state
%%% NW modeling the concentrations we use in our experiments.
%%% NW adding ECS fraction =1-f and osmolytes stuck in the ECS.
%%%%%%% so is the small osmolyte (or impermeable ion) concentration
%%%%%%% zso is its charge
%%%%%%% xo is the large osmolyte (stuck in ECS) concentration
%%%%%%% zo is its charge
%%%%%%% in this case, Cl concentrations in the bath can be different from
%%%%%%% the ECS, because xo is charged. 
%%%%%%
%%%%%%
function [Wnorm,Vmnorm,Wouab,Vmouab]=PLEfn(nao,ko,clo,xo)

%% modeling the addition of sucrose under normal and ouabain-treated conditions



%ouab=[0,0.1,0.5,1]; %%%simulating differnet concentrations of ouabain and their time to effect
ouab=0;

R=25.69*1e-3; F=96485;  % R (RT/F) in Volts, Faraday's constant C/mol
n=1200; % # points to plot 
Vm=1:n-1; K=1:n-1; Na=1:n-1; Cl=1:n-1; W=1:n-1; X=1:n-1; time=1:n-1;%create plotting arrays
time(1)=0; dt=1e-3;  %zero time, dt time step
tt=6000*20; ts=tt/n; tit=round(tt/dt);% tt = total time (s)
ton=0; toff=40000; %time when pump turned on & off 
gna=0.001/F ; gk=0.03/F ; gcl=0.02/F; % conductances in S/dm^2, where dm is decimeters (1/10th of a meter)
%p=0.5e-4/F; %Pump rate (C/dm^2 s) div by F 
p=0.18e-4/F; %Pump rate (C/dm^2 s) div by F 
ck=2;cna=3;%pump stoichiometry
rad=.5; rad=rad*1e-5; %(radius in um convert to dm)
w=(4/3)*pi*rad^3; %vol in liters. spherical cell  
Ar=4*pi*rad^2; %area in dm^2
C=1e-4; % unit membrane capacitance F/dm^2
FinvCAr=F/(C*Ar); %(F/C*Area)

Oso=nao+ko+clo+xo;

x=50e-3; z=-1; MX=x*w; % x is intracellular concentration. MX number of moles of X
cl=(Oso-2*x)/2;   na=0.8*(Oso-cl-x); k=0.2*(Oso-cl-x); %INITIAL INTRA CONCs
ctr=1;t=0; % ctr counter for plotting points, t real time
sw=0;% switch to turn pump on or off
V =FinvCAr*w*(na+k-cl+z*x); % starting voltage

 for i = 2:tit 
     if (toff>t)&&(t>ton) sw=1 ;
  else sw=0;
      if t>toff sw=1*ouab;      
      end 
     end
  V =FinvCAr*w*(na+k-cl+z*x);
  invw=1/w;
  dna=-dt*Ar*invw*(gna*(V-R*log(nao/na))+sw*cna*p);%increase in Na during time step dt
  dk=-dt*Ar*invw*(gk*(V-R*log(ko/k))-sw*ck*p);%increase in K during time step dt
  dcl=dt*Ar*invw*(gcl*(V+R*log(clo/cl)));%increase in Cl during time step dt
  na=na+dna; k=k+dk; cl=cl+dcl; % increment concentrations
  Osi=na+k+cl+x; % intracellular osmolarity 
  w2=(w*Osi)/(Oso); % correct intracellular volume 

  na=(na*w)/w2; k=(k*w)/w2; cl=(cl*w)/w2 ; x=(x*w)/w2; % correct concens 
  w=w2;

  if t>=ctr*ts 
      Vm(ctr)=1000*V; K(ctr)=1000*k; Na(ctr)=1000*na; % save variables for plotting 
      Cl(ctr)=1000*cl; W(ctr)=100*(w/5E-14);X(ctr)=1000*x; time(ctr)=t; ctr=ctr+1;
  end
  t=t+dt; % increment time 
 end
% 
%  figure(1)
%  hold on
%  subplot(2,1,1)
%   hold on
% 
% plot (time,K,'--k',time,Na,'-r', time,Cl,'-g',time,X,'-b')
% legend('K (mM)','Na','Cl','X')
% legend('location','northeastoutside');
% subplot(2,1,2);
%   hold on
% 
% %plot(time,W, time, Vm, time, -W./Vm, '--k');
% plot(time,W/wtot*100, time, Wo/wtot*100, time, Vm);
% legend('intra vol %', 'extracellular vol %','Vm (mV)')
% legend('location','northeastoutside');


Wnorm=W(398);
Vmnorm=Vm(398);
Wouab=W(1198);
Vmouab=Vm(1198);

