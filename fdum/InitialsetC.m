% Author: Rong Wang // contact rongwang@fudan.edu.cn //
% Date: 2021.6.16

global  realtime clim0 Tat airco2 carbonbudget20 ESC INT FFlux

%   clim0:  initial climate variables 1970
%   1 air temperature; 2 ocean temperature; 3 atm carbon (Gt C); 
%   4 upper ocean and biosphere carbon (Gt C); 5 deep ocean carbon (Gt C)
%   Tat: global surface temperature for 1880 - 2020 relative to the 1880-1900 average
%   airco2: atmospheric CO2 levels measured from 1959-2020 at Mauna Loa Observatory, Hawaii
%   carbonbudget20: 1 year; 2 fossil emi; 3 LUC emi; 4 Atm growth; 5 Ocean sink; 6 Land sink; 7 Cement carbonation

% atmospheric CO2 levels measured from 1959-2020 at Mauna Loa Observatory, Hawaii
% https://climate.nasa.gov/vital-signs/carbon-dioxide/
airco2mm=load('files\air_co2_1958_2020.txt');
airco2=zeros(62,1);
for i=1:62
    for j=1:12
        airco2(i,1)=airco2(i,1)+airco2mm((i-1)*12+j+10,5)/12;
    end
end

Tat=load('files\air_temperature_1880_2020.txt');
Tat(:,3)=Tat(:,3) - mean(Tat(1:21,3),1);
% global surface temperature relative to the 1880-1900 average
% https://climate.nasa.gov/vital-signs/global-temperature/

% conversion of a time series to the time format of realtime
airco2 = timeconvert( realtime, airco2, 1959, 1 );
Tat = timeconvert( realtime, Tat, 1880, 1 );

%Global carbon cycle
%Atmospheric carbon (GtC) in 1971 = 860 (for 2010-2019) - air growth for 1971-2014
Mat0 = 860 - sum(carbonbudget20(1:44,4),1);
%Land surface soil carbon (GtC), 203 for active carbon pool from Crowther et al., Nature 2016 550 from Friedlingstein, et al., Global Carbon Budget 2020
Mss0 = Mat0 * sum(carbonbudget20(1:44,6),1) / sum(carbonbudget20(1:44,4),1) * 203/(203+550);
%Land deep soil carbon (GtC), Friedlingstein, et al., Global Carbon Budget 2020
Msd0 = 3600;
%Terrestrial biosphere carbon (GtC), 203 from Crowther et al., Nature 2016 550 from Friedlingstein, et al., Global Carbon Budget 2020
Mbp0 = Mat0 * sum(carbonbudget20(1:44,6),1) / sum(carbonbudget20(1:44,4),1) * 550/(203+550);
%Shallow ocean carbon (GtC)
Mocs0 = Mat0 * sum(carbonbudget20(1:44,5),1) / sum(carbonbudget20(1:44,4),1) ;
%Deep ocean carbon (GtC)
Mocd0 = 1750;
%Ocean biota carbon (GtC)
Mocb0 = 3;
%Initial atmospheric temperature above 1900
Tat0 = Tat(1,3);
%Initial temperature of deep oceans
Tlo0 = 0;
%Dynamic variables in the climate model
%Carbon: 1 air, 2 surface soil, 3 deep soil, 4 terrestrial biosphere, 5 shallow ocean, 6 deep ocean, 7 ocean biota 
%Temperature: 8 surface air, 9 ocean temperature, 10 land carbon, 11, ocean sink
clim0 = [ Mat0, Mss0, Msd0, Mbp0, Mocs0, Mocd0, Mocb0, Tat0, Tlo0, 0, 0 ];

%Carbon flux between pools
% air to surface soil GtC/yr Friedlingstein, GCB 2020
Flux12 = 0; 
% air to land biosphere GtC/yr = NPP from Ciais, NSR, 2020
Flux14 = 50.3; 
% air to surface ocean GtC/yr = NPP from Wang, GRL, 2015
Flux15 = 48.5; 
% land biosphere to soil GtC/yr = soil heteo respiration from Ciais, NSR, 2020
Flux42 = 39.1;
% IPCC-AR5-FAQ 6.2, Figure 1
LTsoil=sqrt(10*500);
% IPCC-AR5-FAQ 6.2, Figure 1
LTocean=sqrt(100*2000);
% surface soil to deep soil GtC/yr
Flux23 = clim0(2)/LTsoil;
% surface ocean to deep ocean GtC/yr
Flux56 = clim0(5)/LTocean;
% deep soil to surface soil GtC/yr
Flux32 = Flux23;
% deep ocean to surface ocean GtC/yr
Flux65 = Flux56;

FFlux = [ESC, INT, Flux12, Flux14, Flux15, Flux42, Flux23, Flux56, Flux32, Flux65];
