% Author: Rong Wang // contact rongwang@fudan.edu.cn //
% Date: 2021.2.10

function [ econcn ]  =  countrydata( output_cap, cndata )

global   alpha elas theta2 EFco2 Egreen
%   output_cap 45x33: 1-23 for model; 24-33 for obs

% Economic variables
% 1 EUE; 2 EPE; 3 ENE; 4 backstop price $/tCO2
% 5 abatement cost as a percentage of GDP; 6 climate damage as a percentage of GDP; 7 net output
% 8 fraction of labor allocation to energy; 9 fraction of investment allocation to energy
% 10 total capital (t$); 11 energy capital (t$); 12 energy (PJ); 13 output (t$); 14 energy price $/kWh; 15 omega
% 16 fraction of energy investment allocated to zero-carbon energy
% 17 energy capital carbon-emission-free (t$);
% 18 fraction to abate CO2 emission; 19 carbon price $/tCO2; 20 CO2 emissions Gt CO2; 21 green energy PJ; 22 invest change, 23 labor change

cn_num=size(cndata,1);
econcn=zeros(49,24,cn_num);
for t=1:49
    for cn=1:cn_num
        %Climate Damage
        D0 = 0;
        %Initial capital stock ($ trillion 2010 USD) 135 in DICE-2013R; 137 in DICE-2007; 223 in DICE-2016R
        K0 = cndata(cn,1*49-48+t);
        %Industrial energy (PJ) 582030
        e0 = cndata(cn,2*49-48+t);
        %Initial gross output (trill 2010 USD)
        q0 = cndata(cn,3*49-48+t)/(1-D0);
        %Initial population (millions)
        L0 = cndata(cn,4*49-48+t);
        %Initial share of energy expenditure in GDP
        se0 = output_cap(min(45,t),32) * (e0/q0) / (cndata(1,2*49-48+t)/cndata(1,3*49-48+t));
        %Industrial emissions (Gt CO2 per year)
        IE0 = e0*EFco2(1,1)*(1-Egreen(min(45,t),8)); % 34.91 in DICE2013
        %Initial level of total factor productivity
        A0 = q0 / (K0^alpha) / (L0/1000)^(1-alpha);
        %Energy price $/kWh
        pe0 = output_cap(min(45,t),31);
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
        %fraction of energy investment allocated to net-zero-carbon green energy
        Vg0 = cndata(cn,5*49-48+min(45,t))/cndata(cn,2*49-48+min(45,t));
        %Capital for carbon-emission-free energy
        Kgreen0 = Ke0 * Vg0;
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
        zc0 = e0*Vg0;
        econcn(t,1:24,cn) = [eue0, epe0, ene0, bs0, abate0, D0, qnet0, Le0, Ve0, K0, Ke0, e0, q0, pe0, se0, Vg0, Kgreen0, acta0, pc0, IE0, zc0,0,0, L0];      
    end
end

end
