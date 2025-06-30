%% Williamson, et al, Magnetic Resonance Letters (2025)
%%% predict cell volumes based on osmotic and ionic condtions
%%% based in cell volume equation # 3 in Kay, Frontiers (2017)
%%% added an osmotic component from trapped extracellular impermeants
%%%inputs: extracellular Na, K, Cl and osmolyte (so), Voltage,  
%%% total volume wtot, 
%%% and trapped intracellular and extracellular impermeants x and xo
function w=cellVolFnPfnf(nao,ko,clo,so,V,wtot,xo,x,z)

R=8.314; F=96485; % R [J/mol K], Faraday's constant [C/mol]=[J/mol volt]
T=273+25; % Temperature in Kelvin
Posmol=nao+ko+clo+so; % osmolarity of extracellular ions and osmolytes

myfun= @(w,z,x,Posmol,xo,wtot,Clo,F,V,R,T) (1-z)*x./(Posmol+xo/(wtot-w)-2*Clo*exp(F*V/(R*T)))-w; 


fun = @(w) myfun(w,z,x,Posmol,xo,wtot,clo,F,V,R,T);    % function of x alone
w = fzero(fun,[0,wtot*.9999])

