%% Williamson, et al, Magnetic Resonance Letters (2025)
%%% predict cell volumes based on osmotic and ionic condtions
%%% This is equation # 3 in Kay, Frontiers (2017)
%%%inputs: extracellular Na, K, Cl and osmolyte (so), Voltage,
%%% and trapped intracellular impermeands

function w=cellVolFn(nao,ko,clo,so,V,x,z)

R=8.314; F=96485; % R (RT/F) in Volts, Faraday's constant C/mol
T=273+25; % Temperature in Kelvin
Posmol=nao+ko+clo+so; % osmolarity of extracellular ions and osmolytes

w=(1-z)*x./(Posmol-2*clo*exp(F*V/(R*T)));
