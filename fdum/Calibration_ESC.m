% Author: Rong Wang // contact rongwang@fudan.edu.cn //
% Date: 2021.6.18

function [ output_esc] = Calibration_ESC( FFlux, displays )

global  realtime Tat airco2 carbonbudget20 clim0 S0
%   realtime:  time
%   clim: climate variables
%   carbon: 1 air, 2 surface soil, 3 deep soil, 4 terrestrial biosphere, 5 shallow ocean, 6 deep ocean, 7 ocean biota 
%   temperature: 8 surface air, 9 ocean temperature
%   flux: 10 land sink GtC/yr; 11 ocean sink GtC/yr
%   output_esc: output of economic variables for model calibration
%   Tat: global surface temperature for 1880 - 2020 relative to the 1880-1900 average
%   airco2: atmospheric CO2 levels measured from 1959-2020 at Mauna Loa Observatory, Hawaii
%   carbonbudget20: 1 year; 2 fossil emi; 3 LUC emi; 4 Atm growth; 5 Ocean sink; 6 Land sink; 7 Cement carbonation

GtC_ppm = airco2(1,1)/clim0(1);

s=zeros(45,17);

% observations
s(1:45,14)=airco2(1:45,1); % air CO2 ppm
s(1:45,15)=Tat(1:45,3); % air temperature deg. C above 1900
s(1:45,16)=carbonbudget20(1:45,6); % Land carbon sink GtC
s(1:45,17)=carbonbudget20(1:45,5); % Ocean carbon sink GtC

t=1;
clim1=clim0;
%Calibration of the initial ocean carbon sink
clim1(2)=clim1(2)*0.98;
clim1(4)=clim1(4)*0.98;
clim1(5)=clim1(5)*0.98;
s(1,1:11)=clim1;
s(1,12)=clim1(1)*GtC_ppm; % air CO2 ppm
s(1,13)=carbonbudget20(1,2); % industrial emission GtC

while realtime(t,1)<2015
    IE=carbonbudget20(t,2) * 3.664 ; % industrial emission Gt CO2    
    clim2=climdyn( t, clim1, FFlux, IE );    
    t=t+1;
    s(t,1:11)=clim2(1,1:11);
    s(t,12)=clim2(1)*GtC_ppm; % air CO2 ppm
    s(t,13)=carbonbudget20(t,2); % industrial emission GtC
    clim1=clim2;
end
s(1,10:11) = s(2,10:11);

% subplot(2,2,1); plot(s(:,12)); hold on; plot(s(:,14));
% subplot(2,2,2); plot(s(:,10)); hold on; plot(s(:,16));
% subplot(2,2,3); plot(s(:,11)); hold on; plot(s(:,17));
% subplot(2,2,4); plot(s(:,8)); hold on; plot(s(:,15));

output_esc = s;

S0(1:45,24:34)=s(1:45,1:11);

if displays==1
    display('calibrated equilibrium sensitivity of climate and fluxes: ESC, Flux12, Flux14, Flux15, Flux42, Flux23, Flux56');
    display(FFlux);
end

% % Calibration of warming in 2020 to 1C according to Figure SPM.2 in IPCC-AR6
% tempref=1;
% output_esc(:,8)=output_esc(:,8)-(output_esc(45,8)-tempref);

end


