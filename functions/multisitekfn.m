%% Williamson, et al, Magnetic Resonance Letters (2025)
% estimate k (i.e. AXR) based on multisite exchange model
% with transmembrane exchange between ECS--ICS 
% and geometric exchange between ECS
% component a is ECS, components b and c are ICS.
function k=multisitekfn(Da,Db,Dc,AXRt,AXRg,fa,fb,fc,outS)

%% 
%%if outS=0 then k is a scalar exchange rate constant
%%if outS=1 the k is set equal to S, a cell with all of the fit information 

format LONGG

R1a=1/1000; %1/ms
R1b=1/1000;% 1/ms
R1c=1/1000;% 1/ms
%R2a=20/1000; %1/ms
%R2b=20/1000;% 1/ms
%R2c=20/1000;% 1/ms
SNR=1000000;
nsim=1;

MixingTimelist=[ 0.2,1,2,4,7,10,20,40,80,160,300]; %spinal cord
%MixingTimelist=10.^(linspace(log10(0.2),log10(400))) %using to show single exp. decay
%MixingTimelist=10.^(linspace(log10(0.001),log10(1000))) %using to show single exp. decay
 taulist1=[0.2,0.593];  % new tauMin=0.2 ms, bmin= 0.089 s/mm2  #point near bs=0, then 3 pts at bs=4.5E9
 taulist2=[0.735,0.58];
% taulist1=[0,0.593];  % new tauMin=0.2 ms, bmin= 0.089 s/mm2  #point near bs=0, then 3 pts at bs=4.5E9
% taulist2=[0.739,0.58];
 b1=taulist1.^3*11.1235;
b2=taulist2.^3*11.1235;
ntau=length(taulist1);
nmix=length(MixingTimelist);

S=struct();
for si=1:length(fa)
    S(si).groundTruthfa=fa(si);
    S(si).scanMat=[1,nsim];
end

for si=1:length(S)

    S(si).MixingTimelist=MixingTimelist;

   % for cii=1:length(S(si).scanMat(:,1))
       % S(si).T(cii).Temperature=cond(cii);
        cii=1;
        k=0;
        for ciii=S(si).scanMat(1,1):S(si).scanMat(end)
                set=ciii;
                k=k+1;
            for i=1:nmix
                for j=1:ntau

                       % I(j)=model2SX(b1(j),b2(j),MixingTimelist(i),De,Di,AXR,R1e,R1i,fi(si));
                      %  I(j)=model3SX(b1(j),b2(j),MixingTimelist(i),Da,Db,Dc,AXR,R1a,R1b,R1c,fa,fb,fc); %model3SX %model1ne
                         [I(j),eigM]=model3SX_2(b1(j),b2(j),MixingTimelist(i),Da,Db,Dc,AXRt,AXRg,R1a,R1b,R1c,fa(si),fb(si),fc(si)); %model3SX %model1ne
                       %I(j)=model3R(taulist1(j),taulist2(j),MixingTimelist(i),De,AXR);
                        I(j)=I(j)+randn(1,1)/SNR;

                end       
                S(si).T(cii).I_2(k,i)=I(1);
                S(si).T(cii).I_3(k,i)=I(2);  
            end
            S(si).T(cii).eigM=eigM;
        end
        % S(si).T(cii).SNR=mean(S(si).T(cii).I_2(:))/mean(S(si).T(cii).noise(:));
   % end
end

S=FitExchangeData_repeatMeasure(S,0); % the DEXR method for estimating exchange.
k=S.T.AXR_I3;
%% overwriting k and outputting the fits (see above)
if outS==1
    k=S;
end

%% define functions
%function [III,I_3]=model3SX(b1,b2,tm,Da,Db,Dc,AXR,R1a,R1b,R1c,fa,fb)
function [III,eigM]=model3SX_2(b1,b2,tm,Da,Db,Dc,AXRt,AXRg,R1a,R1b,R1c,fa,fb,fc)

%fc=1-fa-fb;
DD=[Da,Db,Dc]'.*eye(3);

OD1=expm(-b1*DD);
OD2=expm(-b2*DD);
S0=[fa,fb,fc]';  


KK=     [(fb+fc)*AXRt, -fa*AXRt,         -fa*AXRt; ...
        -fb*AXRt,      fa*AXRt+fc*AXRg, -fb*AXRg; ...
        -fc*AXRt,    -fc*AXRg,         fa*AXRt+fb*AXRg];

%%% check detailed mass balance
%Mdet=KK*S0;
%%%
RR=[R1a,R1b,R1c]'.*eye(3); %%this is the same as model 1 unless R1i not equal to R1e
OE=expm(-tm*(KK+RR));
I_3=OD2*OE*OD1*S0;
III=[1,1,1]*OD2*OE*OD1*S0;
eigM=eig(KK); %% see eq. 3 in Dortch (2012)   https://doi.org/10.1002/mrm.24571 
%%signal decays should be related to eigenvalues...
