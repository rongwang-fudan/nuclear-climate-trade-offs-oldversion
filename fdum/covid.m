% Author: Rong Wang // contact rongwang@fudan.edu.cn //
% Date: 2021.7.28
% Calibtration of efficiency change during COVID-19
%    delta_eff: changes in 1 EUE, 2 EPE, and 3 ENE

function [ scov, LF ] = covid( FFlux, L, iec, calrsav, LR, switcher, covidyear, delta_eff )
%   scov: outputs
%   LF: lost function of the simulations

global  realtime inputs Egreen eypecovid S0
%   realtime:  time
%   Tat: the observed atmospheric temperature for 1971-2015
%   calrsav: 1, EUE0; 2, ENE0; 3, saving rate 1971-2003; 4, saving rate 2004-2008; 5, saving rate 2009-2015
%   Egreen: renewable energy data quad Btu/yr: 1 total; 2 coal; 3 natural gas; 4 oil; 5 nuclear+renewable; 6 nuclear; 7 renewable
%   S0: initital variables
%   eypecovid:  change of E, Y and pe after covid-19

scov=zeros(82,37);
scov(5:29,35:37) = eypecovid(1:25,1:3);

% inertia of the adjustment of labor / investment allocation
inertia=2;

% Run the simulation
t=1;
covid=1;
simus=zeros(130,34);
simus(1,1:34) = S0(1,1:34);
while realtime(t,1)<=2024
    % using the simulation results before covid-19
    if realtime(t,1)<=2015
        econ1 = S0(t,1:23);
        clim1 = S0(t,24:34);
        t=t+1;
        simus(t,1:34) = S0(t,1:34);
        continue;
    end
    %
    if t<33
        rsav=calrsav(3);
    elseif t<38
        rsav=calrsav(4);
    else
        rsav=calrsav(5);
    end
    %Fraction of investment allocated to carbon-emission-free energy
    fracinv=Egreen(min(size(Egreen,1),t+1),8);
    %Energy cost share (Omega) in the past 20 years
    if t<=10
        omega = inputs(1,6)+(inputs(2,6)-inputs(1,6))*(t-11);
    else
        omega = inputs(min(45,t-10),6);
    end
    %Investment for the previous calender year
    investment=0; tt2=0;
    nextyear=floor(realtime(t+1,1)+realtime(t+1,2)/2);
    for i=1:t
        if nextyear==(floor(realtime(i,1)+realtime(i,2)/2)+1)
            investment = investment + rsav * (1 - simus(i,6) - simus(i,5)) * simus(i,13) * realtime(i,2);
            tt2 = tt2 + realtime(i,2);
        end
    end
    investment=investment/tt2;
    %Change of efficiencies in COVID-19
    if realtime(t,1)>covidyear && covid<=11
        deff=[1+delta_eff(1,covid),1+delta_eff(2,covid),1+delta_eff(3,covid)];
        covid=covid+1;
    else
        deff=[1,1,1];
    end
    %Evolution of state from time t to time (t+1)
    D=0;
    econ2 = econdyn(t, L(t), econ1, fracinv, iec, omega, investment, LR, D, inertia, deff, switcher);
    clim2 = climdyn(t, clim1, FFlux, econ2(20));
    t=t+1;
    econ1 = econ2;
    clim1 = clim2;
    simus(t,1:23)=econ2(1,1:23);
    simus(t,24:34)=clim2(1,1:11);
end
scov(1:82,1:34) = simus(45:126,1:34);

LF=zeros(11,4);
LF(:,1) = abs(scov(18:28,12)./mean(scov(5:16,12),1)-scov(18:28,35)); % E
LF(:,2) = abs(scov(18:28,13)./mean(scov(5:16,13),1)-scov(18:28,36)); % Y
LF(:,3) = abs(scov(18:28,14)./mean(scov(5:16,14),1)-scov(18:28,37)); % pe

LF(:,4)=LF(:,1).*LF(:,1)+LF(:,2).*LF(:,2)+LF(:,3).*LF(:,3); % lost function for capital, energy and output

end
