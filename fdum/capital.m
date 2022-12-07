% Author: Rong Wang // contact rongwang@fudan.edu.cn //
% Date: 2021.6.18
% Calibtration of calrsav by capital, energy and ouput
%    savings rate 0.25 for default

function [ scov, LF ] = capital( L, iec, LR, switcher, calrsav )
%   iec:  regression coeffients of efficiencies against omega
%     1-2 for slope/offset; 1-2 for eue/epe - omega; 3 for ene - (1-omega); 4 for ene - omega
%   Damage power coefficient on temperature, 2 in DICE-2016R
%     dpower = 2;
%   Damage coefficient on temperature, 0.00267 in DICE-2013R or 0.00236 in DICE-2016R
%     dcoef = 0.00267;
%   LR: learning rate

global realtime alpha elas EFco2 inputs Egreen econo0
%   realtime:  time
%   Egreen: 1 total; 2 coal; 3 natural gas; 4 oil; 5 nuclear+renewable; 6 nuclear; 7 renewable; 8 fraction of renewable energy
%   Tat: the observed atmospheric temperature for 1971-2015
%   scov: output of economic variables for model calibration
%   calrsav: 1, EUE0; 2, ENE0; 3, saving rate 1971-2003; 4, saving rate 2004-2008; 5, saving rate 2009-2015
%   EFco2: CO2 emission factors for fossil fuel only tCO2 / MJ
%   inputs 45x6: 1 energy PWh; 2 capital trill $; 3 GDP trill $; 4 population mill; 5 energy price ($/kWh); 6 omega
%   econo:  economic variables over time
%     1 EUE; 2 EPE; 3 ENE; 4 backstop price $/tCO2
%     5 abatement cost as a percentage of GDP; 6 climate damage as a percentage of GDP; 7 net output
%     8 fraction of labor allocation to energy; 9 fraction of investment allocation to energy
%     10 total capital t$; 11 energy capital (trill $); 12 energy (PJ); 13 output (t$); 14 energy price $/kWh; 15 omega
%     16 fraction of energy investment allocated to carbon-emission-free energy
%     17 energy capital carbon-emission-free t$
%     18 fraction to abate CO2 emission; 19 carbon price $/tCO2; 20 CO2 emissions Gt CO2
%     21 green energy PJ; 22 invest change, 23 labor change

% inertia of the adjustment of labor / investment allocation
inertia=2;

scov=zeros(45,33);
%Observation
t=1;
while realtime(t,1)<2016
    E = inputs(t,1)*3600; % energy PJ
    K = inputs(t,2); % capital stock ($ trillion 2010 USD)
    Y = inputs(t,3); % gross output (trill 2010 USD)
    P = inputs(t,4); % population (millions)
    pe = inputs(t,5); % $/kWh
    se = inputs(t,6);
    A = Y / (K^alpha) / (P/1000)^(1-alpha); %Initial level of total factor productivity
    eue = se^(elas/(elas-1)) / (E/Y); %Energy use efficiency $ / KJ
    epe = (A^(elas-1) * se)^(1/(elas-1)) / eue; %Energy production efficiency PJ / (trillion $)^0.3 / (billion cap)^0.7
    ene = (A^(elas-1) * (1-se))^(1/(elas-1)); %Non-energy efficiency (t$)^0.7/(billion cap)^0.7    
    IE = E * EFco2(t,1)*(1-Egreen(t,8)); % industrial emission Gt CO2
    scov(t,24:33) = [eue, epe, ene, K, E, Y, P, pe, se, IE];
    t=t+1;
end

scov(1,1:23)=econo0(1,1:23);
scov(1,1)=scov(1,1)*calrsav(1,1); % Calibration of EUE0
scov(1,3)=scov(1,3)*calrsav(1,2); % Calibration of ENE0
scov(1,15) = 1/(1+(scov(1,3)/scov(1,1)/scov(1,2))^(elas-1)); % omega
scov(1,8) = scov(1,15); % fraction of labor allocation to energy
scov(1,9) = scov(1,15); % fraction of investment allocation to energy
scov(1,11) = scov(1,10) *  1/(1+(scov(1,3)/scov(1,1)/scov(1,2))^(elas-1)); % energy capital (t$)

econ1=scov(1,1:23);
t=1;
while realtime(t,1)<2015
    if t<33
        rsav=calrsav(1,3);
    elseif t<38
        rsav=calrsav(1,4);
    else
        rsav=calrsav(1,5);
    end
    %Fraction of investment allocated to carbon-emission-free energy
    fracinv = Egreen(t+1,8);
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
            investment = investment + rsav * (1 - scov(i,6) - scov(i,5)) * scov(i,13) * realtime(i,2);
            tt2 = tt2 + realtime(i,2);
        end
    end
    investment=investment/tt2;
    %Evolution of state from time t to time (t+1)
    deff=[1,1,1];
    D=0;
    econ2 = econdyn(t, L(t), econ1, fracinv, iec, omega, investment, LR, D, inertia, deff, switcher);
    t=t+1;
    scov(t,1:23) = econ2(1,1:23);
    econ1 = econ2;
end

x1=(scov(:,10)-scov(:,27))./std(scov(:,27),[],1); % K
x2=(scov(:,12)-scov(:,28))./std(scov(:,28),[],1); % E
x3=(scov(:,13)-scov(:,29))./std(scov(:,29),[],1); % Y
LF=zeros(1,4);
LF(1)=sum(x1.*x1, 1);
LF(2)=sum(x2.*x2, 1);
LF(3)=sum(x3.*x3, 1);
LF(4)=sum(x1.*x1+x2.*x2+x3.*x3, 1); % lost function for capital, energy and output

end
