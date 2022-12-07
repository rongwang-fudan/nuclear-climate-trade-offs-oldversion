% Author: Rong Wang // contact rongwang@fudan.edu.cn //
% Date: 2021.7.10

function econo2  =  econdyn( t1, Labor, econo1, fracinv, iec, omega, I, LR, D, inertia, deff, switcher )
%   Finds the next state (t1+1) given the current state and actions at time t1
%   fracinv: Fraction of investment allocated to carbon-emission-free energy
%   omega: energy expenditure, abatement cost and climate change damage in the past 20 years
%   IE:  fossil fuel CO2 emission
%   Tatm:  atm temperature
%   dam1d:  coefficients for damage function
%   iec:  rates of induced efficiency changes
%   1-2 for slope/offset; 1-2 for eue/epe - omega; 3 for ene - (1-omega); 4 for ene - omega

global realtime EFco2 Egreen alpha elas theta2 cesin dk_x dk_e dk_green fudanccm_exo
%   realtime:  time
%   econo:  economic variables over time
%   1 EUE; 2 EPE; 3 ENE; 4 backstop price $/tCO2
%   5 abatement cost as a percentage of GDP; 6 climate damage as a percentage of GDP; 7 net output
%   8 fraction of labor allocation to energy; 9 fraction of investment allocation to energy
%   10 capital (trill $); 11 energy capital (trill $); 12 energy (PJ); 13 output (t$); 14 energy price $/kWh; 15 omega
%   16 fraction of energy investment allocated to carbon-emission-free energy
%   17 energy capital carbon-emission-free (trill $);
%   18 fraction to abate CO2 emission; 19 carbon price $/tCO2; 20 CO2 emissions Gt CO2; 21 green energy PJ; 22 invest change, 23 labor change
%   EFco2:  CO2 emission factors for fossil fuel only tCO2 / MJ
%   Egreen: 1 total; 2 coal; 3 natural gas; 4 oil; 5 nuclear+renewable; 6 nuclear; 7 renewable; 8 fraction of renewable energy
%   alpha:  elasticity of output to capital
%   elas:   elasticity of substitution between energy and non-energy in output
%   L:  population (millions)
%   cesin:  input variables 1ke, 2kx, 3eue, 4epe, 5ene, 6dk_x, 7dk_e

%Time step (year)
tstep = realtime(t1,2);

%Economic variables
econo2 = zeros(1,23);

%EUE $ / KJ
if fudanccm_exo==1
    if realtime(t1,1)<2020
        econo2(1) = econo1(1)*(1+iec(17)+iec(2)*(log10(omega)-iec(16)))^tstep;
    else
        % FudanCCM-Endo
        if switcher(4)==1
            econo2(1) = econo1(1)*(1+iec(17)+iec(2)*(log10((omega+min(0.2,econo1(5)))/(1+econo1(6)))-iec(16)))^tstep;
        else
            econo2(1) = econo1(1)*(1+0.004*(1-0.001)^(realtime(t1,1)-2019))^tstep;
        end
    end
else
    econo2(1) = econo1(1)*(1-0.01*(1-0.001)^(realtime(t1,1)-2019))^tstep; % FudanCCM-Exo or DICE
end

%EPE PJ / (trillion $)^0.3 / (billion cap)^0.7
if realtime(t1,1)<2020
    econo2(2) = econo1(2)*(1+iec(21)+iec(7)*(log10(omega)-iec(20)))^tstep;
else
    econo2(2) = econo1(2);
end

%ENE (trillion $)^0.7 / (billion cap)^0.7
if fudanccm_exo==1
    if realtime(t1,1)<2020
        if log10(1-omega)<iec(11)
            econo2(3) = econo1(3)*(1+iec(25)+iec(12)*(log10(1-omega)-iec(24)))^tstep;
        else
            econo2(3) = econo1(3)*(1+iec(27)+iec(14)*(log10(1-omega)-iec(26)))^tstep;
        end
    else
        % FudanCCM-Endo
        if switcher(2)==1
            enex = log10(max(0.6,min(0.99,(1-omega-min(0.2,econo1(5)))))/(1+econo1(6)));
            if enex<iec(11)
                econo2(3) = econo1(3)*(1+iec(25)+iec(12)*(enex-iec(24)))^tstep;
            else
                econo2(3) = econo1(3)*(1+iec(27)+iec(14)*(enex-iec(26)))^tstep;
            end
        else
            econo2(3) = econo1(3)/(1-0.016*exp(-0.002*(realtime(t1,1)-2019)))^tstep;
        end
    end
elseif fudanccm_exo==2
    % FudanCCM-Exo
    econo2(3) = econo1(3)/(1-0.016*exp(-0.001*(realtime(t1,1)-2019)))^tstep;
elseif fudanccm_exo==3
    % DICE
    econo2(3) = econo1(3); % As TFP in DICE
end

%Abate cost as a percentage to output at time t1+1
if fudanccm_exo==3
    econo2(5) = econo1(4) * econo2(1) /theta2/1000 * econo1(18) ^ theta2;
else
    if switcher(8)==1
        econo2(5) = min(0.2, econo1(4) / theta2 * econo1(18)^theta2  * econo1(12) / econo1(13) * EFco2(min(45,t1)) / 1000);
    else
        econo2(5) = 0;
    end
end

%Climate change damage as a percentage to output at time t1+1
if switcher(1)==1
    econo2(6) = D;
else
    econo2(6) = 0;
end

%Economic net output minus damage and abatement
econo2(7) = min(0.999,max(0.001,(1-econo2(6)) * (1-econo2(5)))) * econo1(13);

%Equivalent energy price for a carbon price $/(kJ total energy)
if switcher(3)==1
    if switcher(6)==1
        carbonprice = EFco2(min(45,t1)) * econo1(1,19) / 1000  * (1 - econo1(17) / econo1(11));
    else
        carbonprice = EFco2(min(45,t1)) * econo1(1,19) / 1000;
    end
else
    carbonprice = 0;
end

%Parameters for labor/investment allocation: 1ke, 2kx, 3eue, 4epe, 5ene, 6investment, 7labor, 8timestep, 9carbonprice
cesin = [econo1(11),econo1(10)-econo1(11),econo2(1),econo2(2),econo2(3),I,Labor/1000,tstep, carbonprice ];

%Optimisation of labor/investment allocation
% Method 1 using optimset
% myoptions = optimset('Display','off','FunValCheck','on','algorithm','sqp','MaxFunEvals',100);
% allo = fmincon(@ces_output,[econo1(8) econo1(9)],[],[],[],[],[0.01 0.01],[0.1 0.3], [], myoptions); % return allo // 1 for labor; 2 for investment

% Method 2 direct comparision
steadyomega = 1/(1+(econo2(3)/econo2(1)/econo2(2))^(elas-1)); % omega
if fudanccm_exo==1 || fudanccm_exo==2
    allo = ces_allocation(steadyomega, 1); % 0 for percentile decimal, 1 for thousand decimal
else
    allo = [steadyomega, steadyomega];
end

%Change of labor allocation
if allo(1)>econo1(8)
    econo2(22) = min(allo(1)-econo1(8), max(0.01, econo1(22)*inertia) * tstep ) / tstep;
else
    econo2(22) = max(allo(1)-econo1(8), min(-0.01, econo1(22)*inertia) *tstep ) / tstep;
end

%Allocation of labor to produce energy
if t1==1
    econo2(8) = allo(1);
else
    econo2(8) = econo1(8) + econo2(22) * tstep;
end

%Change of investment allocation
if allo(2)>econo1(9)
    econo2(23) = min(allo(2)-econo1(9), max(0.01, econo1(23)*inertia) * tstep ) / tstep;
else
    econo2(23) = max(allo(2)-econo1(9), min(-0.01, econo1(23)*inertia) *tstep ) / tstep;
end

%Allocation of investment to produce energy
if t1==1
    econo2(9) = allo(2); % initialization
else
    econo2(9) = econo1(9) + econo2(23) * tstep;
end

%Assuming that the allocation reaches steady state immediatly for calculation of SCC in Abatementzone_scc.m
if econo1(23)==-100
    econo2(8) = (econo2(1) * econo2(2))^(elas-1) / ( (econo2(1) * econo2(2))^(elas-1) + econo2(3)^(elas-1) );
    econo2(9) = econo2(8);
end

%Non-energy Capital trill $
Kx = tstep * I * (1-econo2(9)) + (1 - dk_x) ^ tstep * (econo1(10)-econo1(11));

%Energy capital trill $
econo2(11) = tstep * I * econo2(9) + (1 - dk_e) ^ tstep * econo1(11);

%Capital trill $
econo2(10) = Kx + econo2(11);

%Labor for Energy Production
Labor_e = Labor * econo2(8);

%Change of efficiencies due to the pandemic
if fudanccm_exo==1
    econo2(1,1:3) = econo2(1,1:3).*deff;
end

%Energy PJ/yr
if fudanccm_exo==3
    econo2(12) = econo1(12)* (econo2(7) * econo2(1)) / (econo1(7) * econo1(1));
else
    econo2(12) = econo2(2) * econo2(11)^alpha * (Labor_e / 1000)^(1-alpha);
end

%Non-Energy
X = econo2(3) * Kx^alpha * ( (Labor - Labor_e) / 1000)^(1-alpha);

%Output trill$/yr
if fudanccm_exo==3
    econo2(13) = (econo2(10) ^ alpha) * (Labor / 1000)^(1 - alpha) * econo2(3) * deff(3);
else
    econo2(13) = ((econo2(1)*econo2(12))^((elas-1)/elas) + X^((elas-1)/elas))^(1/((elas-1)/elas));
end

%Energy price $/kWh
if fudanccm_exo==1 || fudanccm_exo==2
    econo2(14) = econo2(1) * (1 + (X/econo2(12)/econo2(1))^(1-1/elas))^(1/(elas-1)) *3600;
end

%Share of energy expenditure in GDP
if fudanccm_exo==1 || fudanccm_exo==2
    econo2(15) = econo2(14) / 3600 * econo2(12) / econo2(13);
end

%Fraction of investment allocated to carbon-emission-free energy
econo2(16) = fracinv;

%Capital carbon-emission-free energy (trill$)
econo2(17) = tstep * I * econo2(9) * fracinv + (1-dk_green)^tstep * econo1(17);

%Fraction of emission abatement after excluding Egreen - historical renewable energy without a carbon tax
if fudanccm_exo==3
    econo2(18) = max(0.0001, 1-(1-fracinv)/(1-Egreen(min(45,t1+1),8)));
else
    econo2(18) = max(0.0001, 1-(1-econo2(17)/econo2(11))/(1-Egreen(min(45,t1+1),8)));
end

%Green energy PJ
econo2(21) = econo2(12) * (1-(1-Egreen(min(45,t1+1),8)) * (1-econo2(18)));

%Marginal cost to abate 100% CO2 emissions $/tCO2 based on a learning scheme
econo2(4) = econo1(4);
if realtime(t1,1)>=2010
    econo2(4) = econo2(4) * (min(1+0.05*tstep, max(1e4, econo2(12)) / max(1e4, econo1(12))))^(theta2-1); % Effect of energy expansion on MAC set for 2010 in DICE
end
if realtime(t1,1)>=2020
    if switcher(9)==1
        if switcher(10)==1
            % 20% for Energy<=50PWh; 30% for Energy=100PWh; 40% for Energy>=200PWh
            econo2(4) = econo2(4)*(max(1+0.01*tstep,econo2(21) / econo1(21)))^(log2(1-min(LR*2,max(LR,LR+log (econo1(21)/3600/50) / log(4) * 0.2))));
        else
            econo2(4) = econo2(4)*(max(1+0.01*tstep,econo2(21) / econo1(21)))^(log2(1-LR));
        end
    end
end

%Carbon price $/tCO2
econo2(19) = econo2(4) * econo2(18)^(theta2-1);

%Industrial emission GtCO2/yr
if fudanccm_exo==3
    % use the DICE model 
    econo2(20) = econo2(13) * econo2(1) * (1-econo2(18));
else
    if switcher(7)==1
        econo2(20) = EFco2(min(45,t1)) * econo2(12) * (1-Egreen(min(45,t1+1),8)) * (1-econo2(18));
    else
        econo2(20) = EFco2(min(45,t1)) * econo2(12) * (1-Egreen(min(45,t1+1),8));
    end
end

end



