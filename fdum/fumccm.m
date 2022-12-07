% Author: Rong Wang // contact rongwang@fudan.edu.cn //
% Date: 2021.7.5
tic
clear;
clear global;

global  elas ESC INT FFlux deltarf fudanccm_exo

%Fudan University model for climate change mitigation % 1 fudanccm endo; 2 fudanccm exo; 3 dice
fudanccm_exo=1;
%Elasticity of substitution
elas = 0.4;
%CO2 Forcings of equilibrium CO2 doubling (Wm-2)
deltarf = 3.8;
%Equilibrium sensitivity of climate Sherwood SC, et al.. An Assessment of Earth's Climate Sensitivity Using Multiple Lines of Evidence. Rev Geophys. 2020 Dec;58(4):e2019RG000678.
ESC = 3.1/deltarf;
%Time inertia of climate system to reach equilibirum (year)
INT = 53; % 53 in OSCAR; 38 in DICE
%Climate damage function: dpo for power coefficient on temperature, dcoef for damage as a percentage of GDP for 1 degree warming
dpo = 2;
%Data of climate damage: 1 (high damage, unconfirmed) 1999-Climate change policy_ quantifying uncertainties for damages and optimal carbon taxes; 2 (moderate, used by DICE) 2017-A Survey of Global Impacts of Climate Change: Replication, Survey Methods, and a Statistical Analysis
damagedata = 2;
%Learning rate on the cost curve
LR = 0.2;
%Peak population
LA = 11500;
%Year of COVID-19 outbreak
covidyear = 2020;
%year to initiate mitigation
abtyear = 2025;
%Ultimate fraction of CO2 emission abatement (1 for zero emissions)
abtfrac = 1;
%Time (years) to abate CO2 emissions
abtlen = 10;
% switcher for  C1	C2	S1	S2	S3	S4	S5	T1	T2	T3	T4
switcher = ones(1,10);

%Historical data of economy
Initialset;

%Historical data of climate
InitialsetC;

%Land-use change emissions GtC and radiative forcing for non-CO2
AerosolsLUC;

%Time series of population: L
L = population( LA );

%regression of regional temperature with global temperature in 11 regions
tempz = tempzoneregress( 11, 0 );

%Calibration of climate damage function
[dcoef, xy_damage] = damage( dpo, damagedata );

%Calibration of induced efficiency change
[iec_cn] = Calibration_IEC( cndata, 0 );

%Calibration of equilibrium sensitivity of climate
[output_esc] = Calibration_ESC( FFlux, 1 );

%Calibration of savings rate by capital, energy and ouput
[calrsav, output_cap] = Calibration_CAP( L, iec_cn(1,:), LR, switcher, 1, 0 );

%Calibration of energy price during COVID-19
[oilprice1971_2019, oilprice2020_normalize, energyprice2020, SE, xy_price, co2emi_2019_2021] = energyprice( covidyear );

%Calibration of ENE reduction by COVID-19
[output_covid, deffs] = Calibration_COVID( FFlux, L, iec_cn(1,:), calrsav, LR, switcher, covidyear, 1, 0 );

%Scenarios of mitigation
mcmode=9; % 1 for sens; 2 MC-FUDAM; 3 MC-FUDAM (nonCO2); 4 MC-FUDAM (nonCO2+neg); 5 MC-DICE; 6 MC-DICE (burke); 7 SCC; 8 SCC-nuclear; 9 SCC-nuclear-stochastic; 10 Investment transfers
endsav=1; % 1 for fixed saving rate; 2 for endogeneous saving rate
gdpcn=ModelCN( FFlux, tempz, L, calrsav, output_cap, deffs, covidyear, endsav, mcmode);




