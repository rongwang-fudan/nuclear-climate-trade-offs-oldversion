% Author: Rong Wang // contact rongwang@fudan.edu.cn //
% Date: 2021.6.16

global  realtime inputs alpha elas theta2 econo0 EFco2 Egreen carbonbudget20 eypecovid dk_x dk_e dk_green S0 cou_iform cndata cnenergy

%Capital share
alpha = 0.3;
% rate of non-energy capital depreciation
dk_x=0.1;
% rate of fossil energy capital depreciation
dk_e=0.1;
% rate of green energy capital depreciation
dk_green=0.1;

% time axis from excel
% e.g. 2020 stands for 1 Jan 2020 and 2020.1 stands for 1 Feb 2020
realtime = xlsread('files\time.xlsx');
realtime(:,1) = realtime(:,1) + 0.00001;
T = size(realtime,1);

% 1 id for 222 countries; 2-50 capital stock ($ trillion 2010 USD) from 1971 to 2019; 51-99 total energy (PJ) from 1971 to 2019
% 100-148: gross output (trill 2010 USD) from 1971 to 2019; 149-197 population (millions) from 1971 to 2019; 198-246 green energy (PJ) from 1971 to 2019
load('files\cndata134.dat','-mat');
cndata=cndata2; clear cndata2;

% country information
% 12 regions:
% 1. Total Eastern and Southern Africa; 2. Total Northern Africa; 3. Total Western and Central Africa; 
% 4. Total East Asia; 5. Total South and South-east Asia; 6. Total Western and Central Asia; 7. Total Europe; 8. Total Caribbean; 
% 9. Total Central America; 10. Total North America; 11. Total Oceania; 12. Total South America.
load('files\cou_iform.dat','-mat'); % 1: id for 222 countries; 2: 2 developing/ 1 developed; 3: 12 region id; 4 OECD; 5 id for 112 countries; 6 pi temperature
cou_iform=zeros(222,7);
cou_iform(:,1:6)=cou_iform2(:,1:6); clear cou_iform2;
%1  N Ame; 2 S Ame; 3 W Eur; 4 N Afr & Middle East; 5 S Afr; 6 Russia; 7 E Asia; 8 S & SE Asia; 9 Oceania; 10 Polar; 11 Ocean
for cn=1:222
    rid=cou_iform(cn,3);
    if rid==10
        cou_iform(cn,7)=1;
    elseif rid==8 || rid==9 || rid==12
        cou_iform(cn,7)=2;
    elseif rid==7 && cn~=165
        cou_iform(cn,7)=3;
    elseif rid==2 || rid==6
        cou_iform(cn,7)=4;
    elseif rid==1 || rid==3
        cou_iform(cn,7)=5;
    elseif cn==165
        cou_iform(cn,7)=6;
    elseif rid==4
        cou_iform(cn,7)=7;
    elseif rid==5
        cou_iform(cn,7)=8;
    elseif rid==11
        cou_iform(cn,7)=9;
    end
end

% all variables
S0=zeros(T,32); % 23+11
% 45x31; 1 energy price ($/kWh) 2-7 energy per capita (MWh/cap) 8-13 gdp per capita (k$/cap) globe china us eu japan india 
% 14-19 omega 20-25 capital stock per capita ($/cap) 26-31 population (million)
% load('files\EconomicHistoricalData_1971_2015.dat','-mat');
% inputs=input; clear input;

% inputs 45x6 from IEA: 1 energy PWh; 2 capital trill $; 3 GDP trill $; 4 population mill; 5 energy price ($/kWh); 6 omega
inputs = load('files\global economic data 1971-2015.txt');

% change of E, Y and pe after covid-19
eypecovid = load('files\covid_E_Y_pe.txt');

% total CO2 emissions from fossil fuel and cements (Gt CO2 yr-1)
% IEco2=load('files\co2emi_1750_2019.txt')/1000;
% https://github.com/owid/co2-data/blob/master/owid-co2-codebook.csv

% 1 year; 2 fossil emi GtC; 3 LUC emi GtC; 4 Atm growth GtC; 5 Ocean sink GtC; 6 Land sink GtC; 7 Cement carbonation
% Budget = 2 + 3 - 4 - 5 - 6 - 7
% Friedlingstein, et al., Global Carbon Budget 2020, https://doi.org/10.5194/essd-12-3269-2020.
carbonbudget20 = load('files\global_carbon_budget_1959_2019.txt');

% renewable energy data quad Btu/yr: 1 total; 2 coal; 3 natural gas; 4 oil; 5 nuclear+renewable; 6 nuclear; 7 renewable
% data from eia (available by country) https://www.eia.gov/international/data/world
% Annual Total energy production INT-Export-04-22-2021_16-49-50.csv
genep=load('files\global energy production 1980-2018.txt');

% fraction of renewable energy
Egreen=zeros(45,8);
% CO2 emission factors for fossil fuel only tCO2 / MJ
EFco2=zeros(45,1); % tCO2 / MJ
for i=1:45
    Egreen(i,1:7)=genep(1:7,min(39,max(1,i-9)));
    Egreen(i,8)=Egreen(i,7)/(Egreen(i,1)-Egreen(i,6)); % excluding nuclear
    EFco2(i,1)=carbonbudget20(i+12,2) * 3.664 / (inputs(i,1)*3600*(1-Egreen(i,8))); % tCO2 / MJ
end
clear genep;

% energy production by country
% 1	Production (quad Btu)
% 2	Coal (quad Btu)
% 3	Natural gas (quad Btu)
% 4	Petroleum and other liquids (quad Btu)
% 5	Nuclear, renewables, and other (quad Btu)
% 6	Nuclear (quad Btu)
% 7	Renewables and other (quad Btu)
cnenergy2=load('files\energyproductionbycountry.txt'); % 193 (countries) * 7 (sources) * 42 (years 1980-2019)
cnenergy=zeros(135,49,7); % 1 for globe; 2-135 for countries
for i=1:193
    cn2=cnenergy2((i-1)*7+1,1);
    cn3=find(cndata(:,1)==cn2);
    if size(cn3,1)==0
        continue;
    end
    for j=1:7
        cnenergy(cn3,10:49,j)=cnenergy2((i-1)*7+j,3:42);
        if cnenergy(cn3,49,j)<=0
            cnenergy(cn3,49,j)=cnenergy2(j,42)*cndata(cn3,4*49-48+49)/cndata(1,4*49-48+49); % interpolate by population
        end
        idx3=find(cnenergy(cn3,1:48,j)<=0);
        for t=idx3(1,end):-1:1
            cnenergy(cn3,t,j)=cnenergy(cn3,t+1,j)*cndata(cn3,4*49-48+t)/cndata(cn3,4*49-48+t+1); % interpolate by population
        end
    end
end
clear cnenergy2;

% conversion of a time series to the time format of realtime
inputs = timeconvert( realtime, inputs, 1971, 1 );
carbonbudget20 = timeconvert( realtime, carbonbudget20, 1959, 1 );
Egreen = timeconvert( realtime, Egreen, 1971, 1 );
EFco2 = timeconvert( realtime, EFco2, 1971, 1 );

%Coefficient in abatement cost curve
theta2 = 2.6;
%Industrial energy (PJ) 582030
e0 = inputs(1,1)*3600; % PWh -> PJ
%Initial share of energy expenditure in GDP
se0 = inputs(1,6);
%Industrial emissions (Gt CO2 per year)
IE0 = e0*EFco2(1,1)*(1-Egreen(1,8)); % 34.91 in DICE2013
%Initial capital stock ($ trillion 2010 USD) 135 in DICE-2013R; 137 in DICE-2007; 223 in DICE-2016R
K0 = inputs(1,2);
%Climate Damage
D0 = 0;
%Initial population (millions)
L0 = inputs(1,4);
%Initial gross output (trill 2010 USD)
q0 = inputs(1,3) / (1-D0);
%Initial level of total factor productivity
A0 = q0 / (K0^alpha) / (L0/1000)^(1-alpha);
%Energy price $/kWh
pe0 = inputs(1,5);
%Energy use efficiency $ / KJ
eue0 = se0^(elas/(elas-1)) / (e0/q0); 
%Energy production efficiency PJ / (trillion $)^0.3 / (billion cap)^0.7
epe0 = (A0^(elas-1) * se0)^(1/(elas-1)) / eue0; 
%Non-energy efficiency (trillion $)^0.7 / (billion cap)^0.7
ene0 = (A0^(elas-1) * (1-se0))^(1/(elas-1));
%Capital for Energy Production
Ke0 = K0 * se0;
%Labor allocation to energy
Le0 = se0;
%Investment allocation to energy
Ve0 = se0;
%Capital for carbon-emission-free energy
Kgreen0 = Ke0 * Egreen(1,8);
%Carbon-emission-free price $/tCO2
bs0 = 550;
%Fraction of CO2 abatements
acta0 = 0;
%Carbon price $/tCO2
pc0 = bs0 * acta0^(theta2-1);
%Abatement cost as a percentage of GDP
abate0 = bs0 / theta2 * acta0^theta2  * e0 / q0 * EFco2(1) / 1000;
%Net output
qnet0 = q0 * (1-abate0-D0);
%green energy PJ
egc0 = e0 * Egreen(1,8);

%Dynamic variables in the economic model
% 1 EUE; 2 EPE; 3 ENE; 4 backstop price $/tCO2
% 5 abatement cost as a percentage of GDP; 6 climate damage as a percentage of GDP; 7 net output
% 8 fraction of labor allocation to energy; 9 fraction of investment allocation to energy
% 10 total capital (t$); 11 energy capital (t$); 12 energy (PJ); 13 output (t$); 14 energy price $/kWh; 15 omega
% 16 fraction of energy investment allocated to carbon-emission-free energy
% 17 energy capital carbon-emission-free (t$);
% 18 fraction to abate CO2 emission; 19 carbon price $/tCO2; 20 CO2 emissions Gt CO2; 21 green energy PJ; 22 invest change, 23 labor change
econo0 = [eue0, epe0, ene0, bs0, abate0, D0, qnet0, Le0, Ve0, K0, Ke0, e0, q0, pe0, se0, Egreen(1,8), Kgreen0, acta0, pc0, IE0, egc0,0,0,1,1];
