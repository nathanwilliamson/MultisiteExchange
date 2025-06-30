function S=FitExchangeData_repeatMeasure(S,plot)
%%% branches from diffusion_exchange_rate Supplementary data and code for: Real-time measurement of diffusion exchange rate in biological tissue; Williamson et al. Journal of Magnetic Resonance (2020)

%% Initialization.
%clear S
%clc
%close all hidden
%% Random stream.
random_seed = round( sum( 1e6 * clock() ) );
s = RandStream('mt19937ar', 'Seed', random_seed);
RandStream.setGlobalStream(s);
%load('dataHold.mat');
%%
MixingTimelist=S.MixingTimelist; %spinal cord
%% fit exchange rate
addpath('diffusion_semiparametric');
number_of_fits = 100;
number_of_mc_fits = 0;
model = {{'exponential'}};
baseline = true; 
clear I_AXR
MixingTimepts=[1:length(MixingTimelist)];
for si=1: length(S)
%for cii=1:length(S(si).scanMat(:,1))
cii=1;
k=0;
for ciii=S(si).scanMat(1,1):S(si).scanMat(end)
k=k+1;
%% fit exchange rate with I(2) normalization
S(si).T(cii).I3_N24(k,:)=S(si).T(cii).I_3(k,:)./S(si).T(cii).I_2(k,:);
%I_AXR_I123(k,:)=I_AXR_I123(k,:)/I_AXR_I123(k,1)
%%% Fit model and estimate parameters.
S(si).T(cii).fit_I3_N24{k,1} = analyze(MixingTimelist(MixingTimepts), S(si).T(cii).I3_N24(k,:), model, true, number_of_fits, number_of_mc_fits);
%plot_fit_and_residuals(MixingTimelist(MixingTimepts), S(si).T(cii).I3_N24(k,:), S(si).T(cii).fit_I3_N24{k,1});
S(si).T(cii).I0_I3_N24(k,1)=S(si).T(cii).fit_I3_N24{k,1}.I0;
S(si).T(cii).baseline_I3_N24(k,1)=S(si).T(cii).fit_I3_N24{k,1}.baseline;
S(si).T(cii).AXR_I3_N24(k,1)=S(si).T(cii).fit_I3_N24{k,1}.components{1,1}.D;
%% fit apparent T1 from I(2), 
%%% Fit model and estimate parameters.
S(si).T(cii).fit_I2{k,1} = analyze(MixingTimelist(MixingTimepts), S(si).T(cii).I_2(k,:), model, false, number_of_fits, number_of_mc_fits);
if plot==1
% plot_fit_and_residuals(MixingTimelist(MixingTimepts), S(si).T(cii).I_2(k,:), S(si).T(cii).fit_I2{k,1});
end
S(si).T(cii).I0_T1app2(k,1)=S(si).T(cii).fit_I2{k,1}.I0;
S(si).T(cii).R1_T1app2(k,1)=S(si).T(cii).fit_I2{k,1}.components{1,1}.D;
% 
%% fit exchange rate from I(3), accounting for T1 in the tissue
S(si).T(cii).I_AXR_I3(k,:)= S(si).T(cii).I_3(k,MixingTimepts)./exp(-MixingTimelist(MixingTimepts)*S(si).T(cii).R1_T1app2(k,1));
%%% HERE changed 20221014% S(si).T(cii).I_AXR_I3(k,:)=I_AXR_I3_temp(k,MixingTimepts)-mean(I_AXR_I3_temp(k,MixingTimeptsBaseline),2);
%%% S(si).T(cii).I_AXR_I3(k,:)=I_AXR_I3_temp(k,MixingTimepts);
%%% Fit model and estimate parameters.
S(si).T(cii).fit_I3{k,1} = analyze(MixingTimelist(MixingTimepts), S(si).T(cii).I_AXR_I3(k,:), model, true, number_of_fits, number_of_mc_fits);
if plot==1
plot_fit_and_residuals(MixingTimelist(MixingTimepts), S(si).T(cii).I_AXR_I3(k,:), S(si).T(cii).fit_I3{k,1} );
end 
S(si).T(cii).I0_I3(k,1)=S(si).T(cii).fit_I3{k,1}.I0;
S(si).T(cii).AXR_I3(k,1)=S(si).T(cii).fit_I3{k,1}.components{1,1}.D;
S(si).T(cii).baseline_I3(k,1)=S(si).T(cii).fit_I3{k,1}.baseline;
%% fit apparent T1 from I(2), with biexp model 
%%% Fit model and estimate parameters.
S(si).T(cii).fit_I2biexp{k,1} = analyze(MixingTimelist(MixingTimepts), S(si).T(cii).I_2(k,:), {{'exponential'},{'exponential'}}, false, number_of_fits, number_of_mc_fits);
if plot==1
plot_fit_and_residuals(MixingTimelist(MixingTimepts), S(si).T(cii).I_2(k,:), S(si).T(cii).fit_I2biexp{k,1});
end
S(si).T(cii).I0_T1app2biexp(k,1)=S(si).T(cii).fit_I2biexp{k,1}.I0;
if S(si).T(cii).fit_I2biexp{k,1}.components{1,1}.theta>S(si).T(cii).fit_I2biexp{k,1}.components{2,1}.theta
S(si).T(cii).R1_T1app2biexp_c1(k,1)=S(si).T(cii).fit_I2biexp{k,1}.components{1,1}.D;
S(si).T(cii).R1_T1app2biexp_c2(k,1)=S(si).T(cii).fit_I2biexp{k,1}.components{2,1}.D;
S(si).T(cii).R1_T1app2biexp_theta(k,1)=S(si).T(cii).fit_I2biexp{k,1}.components{1,1}.theta;
else
S(si).T(cii).R1_T1app2biexp_c1(k,1)=S(si).T(cii).fit_I2biexp{k,1}.components{2,1}.D;
S(si).T(cii).R1_T1app2biexp_c2(k,1)=S(si).T(cii).fit_I2biexp{k,1}.components{1,1}.D;
S(si).T(cii).R1_T1app2biexp_theta(k,1)=S(si).T(cii).fit_I2biexp{k,1}.components{2,1}.theta;
end
S(si).T(cii).R1_T1app2biexp_mean=S(si).T(cii).R1_T1app2biexp_c1(k,1)*S(si).T(cii).R1_T1app2biexp_theta(k,1) + S(si).T(cii).R1_T1app2biexp_c2(k,1)*(1-S(si).T(cii).R1_T1app2biexp_theta(k,1));
%% fit exchange rate from I(3), accounting for T1 in the tissue
S(si).T(cii).I_AXR_I3_subbi(k,:)= S(si).T(cii).I_3(k,MixingTimepts)./(S(si).T(cii).R1_T1app2biexp_theta(k,1)*exp(-MixingTimelist(MixingTimepts)*S(si).T(cii).R1_T1app2biexp_c1(k,1))+(1-S(si).T(cii).R1_T1app2biexp_theta(k,1))*exp(-MixingTimelist(MixingTimepts)*S(si).T(cii).R1_T1app2biexp_c2(k,1)));
S(si).T(cii).fit_I3_subbi{k,1} = analyze(MixingTimelist(MixingTimepts), S(si).T(cii).I_AXR_I3_subbi(k,:), model, true, number_of_fits, number_of_mc_fits);
if plot==1
plot_fit_and_residuals(MixingTimelist(MixingTimepts), S(si).T(cii).I_AXR_I3_subbi(k,:), S(si).T(cii).fit_I3_subbi{k,1} );
end
S(si).T(cii).I0_I3_subbi(k,1)=S(si).T(cii).fit_I3_subbi{k,1}.I0;
S(si).T(cii).AXR_I3_subbi(k,1)=S(si).T(cii).fit_I3_subbi{k,1}.components{1,1}.D;
S(si).T(cii).baseline_I3_subbi(k,1)=S(si).T(cii).fit_I3_subbi{k,1}.baseline;
%% fit exchange rate from I(3), accounting for T1 in the tissue, with biexp model
% 
% S(si).T(cii).fit_I3_bibi{k,1} = analyze(MixingTimelist(MixingTimepts), S(si).T(cii).I_AXR_I3_subbi(k,:), {{'exponential'},{'exponential'}}, true, number_of_fits, number_of_mc_fits);
% 
% 
% if plot==1
% plot_fit_and_residuals(MixingTimelist(MixingTimepts), S(si).T(cii).I_AXR_I3_subbi(k,:), S(si).T(cii).fit_I3_subbi{k,1} );
% end
% S(si).T(cii).I0_I3_bibi(k)=S(si).T(cii).fit_I3_bibi{k,1}.I0;
% if S(si).T(cii).fit_I3_bibi{k,1}.components{1,1}.theta>S(si).T(cii).fit_I3_bibi{k,1}.components{2,1}.theta
% S(si).T(cii).AXR_I3_bibi_c1(k,1)=S(si).T(cii).fit_I3_bibi{k,1}.components{1,1}.D;
% S(si).T(cii).AXR_I3_bibi_c2(k,1)=S(si).T(cii).fit_I3_bibi{k,1}.components{2,1}.D;
% S(si).T(cii).AXR_I3_bibi_theta(k,1)=S(si).T(cii).fit_I3_bibi{k,1}.components{1,1}.theta;
% else
% S(si).T(cii).AXR_I3_bibi_c1(k,1)=S(si).T(cii).fit_I3_bibi{k,1}.components{2,1}.D;
% S(si).T(cii).AXR_I3_bibi_c2(k,1)=S(si).T(cii).fit_I3_bibi{k,1}.components{1,1}.D;
% S(si).T(cii).AXR_I3_bibi_theta(k,1)=S(si).T(cii).fit_I3_bibi{k,1}.components{2,1}.theta;
% end
% S(si).T(cii).AXR_I3_bibi_mean=S(si).T(cii).AXR_I3_bibi_c1(k,1)*S(si).T(cii).AXR_I3_bibi_theta(k,1) + S(si).T(cii).AXR_I3_bibi_c2(k,1)*(1-S(si).T(cii).AXR_I3_bibi_theta(k,1));
% 
end
end
%save('Data_SpinalCord.mat','S');
