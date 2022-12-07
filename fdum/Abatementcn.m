% Author: Rong Wang // contact rongwang@fudan.edu.cn //
% Date: 2021.7.31
% Establishing the energy system by zone

function [ npv, tip, S ] = Abatementcn( FFlux, L, dpower, dcoef, tempopt, LR, covidyear, iec_cn, calrsav, deffs, techdiff, switcher, econcn, abtcn, tempdiff, weitzman, dicetfp, tipping, emi_neg, optlrsav, elasmu, rour, scc )

global Egreen realtime S0 cou_iform theta2 alpha fudanccm_exo EFco2 burkecoef tempzone cndata cnenergy nuclear_tip nuclear_tech nuclear_acc randnumc invest_trans energy_growth mac_curves neg_year

% cou_iform % 1: id for 222 countries; 2: 2 developing/ 1 developed; 3: region id; 4 OECD; 5 id for 112 countries; 6 pi temperature

T = size(realtime,1);

% inertia of the adjustment of labor / investment allocation
inertia=2;

% initialize the DICE model
cn_num=size(cndata,1);
if fudanccm_exo==3
    % DICE
    TFPd = zeros(cn_num,49); % Total factor productivity in the DICE model
    EMId = zeros(cn_num,49); % Emission Intensity (emission/Y)
    for t=1:49
        for cn=1:cn_num
            TFPd(cn,t) = econcn(t,13,cn)/(econcn(t,10,cn)^alpha)/(cndata(cn,4*49-48+t)/1000)^(1-alpha);
            EMId(cn,t) = econcn(t,20,cn)/econcn(t,13,cn);
        end
    end
end

abatecn2019=zeros(cn_num,1);
heatratio = sum(cnenergy(1,49,2:4),3) / cnenergy(1,49,5) * 0.28 / (1-0.28);
for cn=1:cn_num
    if cnenergy(cn,49,1)>0
        abatecn2019(cn,1)=cnenergy(cn,49,5)./(sum(cnenergy(cn,49,2:4),3)/heatratio+cnenergy(cn,49,5));
    else
        abatecn2019(cn,1)=cnenergy(1,49,5)./(sum(cnenergy(1,49,2:4),3)/heatratio+cnenergy(1,49,5));
    end
end

% 36 utility; 37-44 for nuclear: 37 power GW; 38 PWh; 39 CO2 emissions; 40 capacity factor; 41 number of accident; 42 damage as a % to GDP; 43 nuclear costs t$; 44 number of accident level>=5
% 45 climate investment contributed by a country; 46 climate investment used in a country; 47 emission abatements Gt CO2; 48 national carbon price $/tCO2
% 49-60 climate investment transfer from region 1-12
if invest_trans(3)>0
    endyear=2201;
    S=zeros(T,60,cn_num+1); 
else
    endyear=2301;
    S=zeros(T,44,cn_num+1); 
end

S(1,1:34,1)=S0(1,1:34);
S(1,35,1)=L(1);
S(1,1:23,2:(cn_num+1))=econcn(1,1:23,1:cn_num);
S(1,35,2:(cn_num+1))=cndata(1:cn_num,4*49-48+1); % Labor in 1971

% Tipping Interaction according to Cai 2016 NCC: 1 for Atlantic meridional overturning circulation (AMOC), 2 for melt of the Greenland ice sheet (GIS), 
% 3 for collapse of the West Antarctic ice sheet (WAIS), 4 for retreat of the Amazon rainforest (AMAZ), and 5 for shift to a more persistent El Niño regime (ENSO)
tip_para=[0.063 50 15 0; 0.188 1500 10 0.067; 0.104 500 5 0.2; 0.163 50 5 1; 0.053 50 10 0.2];
tip_fij=[0 -0.235 0.125 0.55 0.121; 1.62 0 0.378 0.108 0; 0.107 0.246 0 0 0; 0 0 0 0 0; -0.083 0 0.5 2.059 0];
tip_prob=ones(5,4);
tip=ones(T,11);
t=1;
covid=1;
ZCE0 = zeros(5,2);
id_nuclear = zeros(cn_num-1,2);
idddd=0;

while realtime(t,1)<endyear
    % using the simulation results before abatements
    if realtime(t,1)<2018 % using the historical data
        econ1 = S0(t,1:23);
        clim1 = S0(t,24:34);
        t=t+1;
        S(t,1:34,1) = S0(t,1:34);
        S(t,35,1) = L(t);
        S(t,1:23,2:(cn_num+1)) = econcn(t,1:23,1:cn_num);
        S(t,4,2:(cn_num+1))  = S0(t,4); % MAC of zero-carbon energy 
        S(t,22,2:(cn_num+1))  = S0(t,22);
        S(t,23,2:(cn_num+1))  = S0(t,23);
        S(t,35,2:(cn_num+1)) = econcn(t,24,1:cn_num);        
        %Initializing TFP adn Emission_factor in DICE
        if fudanccm_exo==3
            S(t,1,1)=EMId(1,floor(realtime(t,1))-1970);
            S(t,3,1)=TFPd(1,floor(realtime(t,1))-1970);
            S(t,1,2:(cn_num+1))=EMId(1:cn_num,floor(realtime(t,1))-1970);
            S(t,3,2:(cn_num+1))=TFPd(1:cn_num,floor(realtime(t,1))-1970);
        end
        continue;
    end
    
    %Endogeneous saving rate
    if t<33
        rsav=calrsav(3);
    elseif t<38
        rsav=calrsav(4);
    elseif t<121
        rsav=calrsav(5);
    else
        rsav=optlrsav(t-120);
    end
    
    %Fraction of investment allocated to carbon-emission-free energy: S transition
    cn_yabt=abtcn(1,1);
    cn_fabt=abtcn(1,2);
    if floor(realtime(t,1))>=cn_yabt
        cn_alen=abtcn(1,3);
    else
        cn_alen=100;
    end
    cn_t0=cn_yabt-cn_alen/100*(cn_yabt-2025);
    cn_Nabt=max(0,realtime(t,1)-cn_t0);
    fracinv=(Egreen(45,8)-cn_fabt)*exp(-(cn_Nabt^2)/2/cn_alen/cn_alen)+cn_fabt; 
    
    %Energy cost share (Omega) in the past 20 years
    omega=0; tt=0;
    for t2=1:t
        if (realtime(t,1)-realtime(t2,1))<20
            omega=omega+S(t2,15,1)*realtime(t2,2);
            tt=tt+realtime(t2,2);
        end
    end
    omega = omega/tt;
    
    %Investment for the previous calender year
    investment=0; tt2=0;
    nextyear=floor(realtime(t+1,1)+realtime(t+1,2)/2);
    for i=1:t
        if nextyear==(floor(realtime(i,1)+realtime(i,2)/2)+1)
            investment = investment + rsav * S(i,7,1) * realtime(i,2);
            tt2 = tt2 + realtime(i,2);
        end
    end
    investment=investment/tt2;
    
    %Change of efficiencies in COVID-19
    if realtime(t,1)>covidyear && covid<=11 && realtime(t,1)<(covidyear+1)
        deff=[1+deffs(1,covid),1+deffs(2,covid),1+deffs(3,covid)];
        covid=covid+1;
    else
        deff=[1,1,1];
    end
    
    %Global warming
    tempw = clim1(8); % warming in year t relative to 1850-1990
    tempw0 = S0(49,31); % warming in 2019 relative to 1850-1990
    
    %Tipping function according to Cai 2016 NCC
    for tip_i=1:5
        %Interaction
        interact=0;
        for tip_j=1:5
            if tip_prob(tip_j,2)==0 && tipping==2
                interact=interact+tip_fij(tip_j,tip_i);
            end
        end
        if nuclear_tip==0
            tip_prob(tip_i,2)=tip_prob(tip_i,2)*exp(-tip_para(tip_i,1)/100*max(0,tempw-1)*(1+interact)); % Probabiity of no tipping
            if tip_prob(tip_i,2)<0.9
                tip_prob(tip_i,2)=0; % Record the time that a tipping might occur when the probablity is >10%
                if tip_prob(tip_i,4)==1
                    tip_prob(tip_i,4)=t; % Time to reach a tipping point
                end
            end
            tip(t,tip_i)=tip(t-1,tip_i)*exp(-tip_para(tip_i,1)/100*max(0,tempw-1)*(1+interact));
        elseif nuclear_tip==1 || nuclear_tip==2
            if tempw>=2
                tip_prob(tip_i,2)=tip_prob(tip_i,2)*exp(-tip_para(tip_i,1)/100*(tempw-1)*2*(1+interact)); % Probabiity of no tipping
            else
                tip_prob(tip_i,2)=tip_prob(tip_i,2)*exp(-tip_para(tip_i,1)/100*max(0,tempw-1)*(1+interact)); % Probabiity of no tipping
            end
            if tip_prob(tip_i,2)<randnumc(tip_i) % uniform distribution
                tip_prob(tip_i,2)=0; % Record the time that the probablity the tipping exceeds a random value in uniform distribution 0-1
                if tip_prob(tip_i,4)==1
                    tip_prob(tip_i,4)=t; % Time to reach a tipping point
                end
                tip(t,tip_i)=0; % tipping occurs
            else
                if tempw>=2
                    tip(t,tip_i)=tip(t-1,tip_i)*exp(-tip_para(tip_i,1)/100*(tempw-1)*2*(1+interact));
                else
                    tip(t,tip_i)=tip(t-1,tip_i)*exp(-tip_para(tip_i,1)/100*max(0,tempw-1)*(1+interact));
                end
            end
        end
    end
    
    %Climate change catastrophe according to Weitzman 2012
    if weitzman==1
        D = 1 - (1+dcoef * tempw0^dpower+(tempw0/6.081)^6.754)/(1+dcoef * tempw^dpower+(tempw/6.081)^6.754);
    else
        D = 1 - (1+dcoef* tempw0^dpower)/(1+dcoef * tempw^dpower);
    end
    
    %Tipping function according to Cai 2016 NCC
    damage_tip=1;
    if tipping>0
        for tip_i=1:5
            if tip_prob(tip_i,2)==0
                if nuclear_tip==0
                    damage_tip_i = tip_para(tip_i,3)/100 * min( 1,(t-tip_prob(tip_i,4))/tip_para(tip_i,2) ) ;
                elseif nuclear_tip==1 || nuclear_tip==2
                    damage_tip_i = tip_para(tip_i,3)/100;
                end
                damage_tip = damage_tip * (1-damage_tip_i * (1-tip(t,tip_i)));
                tip(t,tip_i+5) = damage_tip_i * (1-tip(t,tip_i));
            end
        end
    end
    tip(t,11) = 1-damage_tip;
    D = 1-(1-D)*damage_tip;
    
    %Economy
    econ2 = econdyn(t, L(t), econ1, fracinv, iec_cn(1,:), omega, investment, LR, D, inertia, deff, switcher);
    
    %National loop
    for cn=1:cn_num
        cn_labor = S(t,35,cn+1);
        cn_econ0 = S(t,1:23,cn+1);
        cn_omega = S(t,15,cn+1);
        cn_invest = rsav*S(t,7,cn+1);
        cn_growth=abtcn(cn,4); % the prescribed initital rate of per capita GDP growth in 2025
        
        % Investment allocation
        cn_yabt=abtcn(cn,1);
        cn_fabt=abtcn(cn,2);
        if floor(realtime(t,1))>=cn_yabt
            cn_alen=abtcn(cn,3);
        else
            cn_alen=100;
        end
        cn_t0=cn_yabt-cn_alen/100*(cn_yabt-2025);
        cn_Nabt=max(0,realtime(t,1)-cn_t0);
        cn_fracinv=(cndata(cn,5*49-48+45)/cndata(cn,2*49-48+45)-cn_fabt)*exp(-(cn_Nabt^2)/2/cn_alen/cn_alen)+cn_fabt;
        
%         cn_fracinv= max(0, min(1,cn_climate_invest / (cn_invest * S(t,9,cn+1))));
        
        % Absolute temperature
        cntemp0=tempopt; % global temperature in 2016-2019
        cntempw0=tempw0; % global warming in 2019 relative to 1850-1990
        cntempw=tempw;   % global warming in year t relative to 1850-1990
        if tempdiff==3 && cn>1
            cnzone=cou_iform(cndata(cn,1),7);
            cntempw0=tempw0*tempzone(cnzone,1)+tempzone(cnzone,2); % warming by country in 2016-2019
            cntempw=tempw*tempzone(cnzone,1)+tempzone(cnzone,2); % warming by country in year t
        end
        if tempdiff>1 && cn>1
            cntemp0=cou_iform(cndata(cn,1),6); % absolute temperature by country in 2016-2019
        end
        cntemp=cntemp0-cntempw0+cntempw; % absolute temperature by country in year t
        
        % Damage
        if weitzman==1
            D = 1 -(1+dcoef*(abs(cntemp0-tempopt))^dpower+(tempw0/6.081)^6.754) /(1+dcoef*(abs(cntemp-tempopt))^dpower+(tempw/6.081)^6.754);
        else
            D = 1 -(1+dcoef*(abs(cntemp0-tempopt))^dpower) /(1+dcoef*(abs(cntemp-tempopt))^dpower);
        end
        D = 1 -(1-D)*damage_tip;
        
        % Calibration of iec to match the growth of economy in the future
        if t<=121 && cn>1
            % EUE
            iec_cn(cn,16)=log10(cn_omega);
            iec_cn(cn,17)=0.01; % like DICE
            iec_cn(cn,28)=cn_omega;
            % ENE
            enerate=((1-0.3)*cn_growth-cn_omega*iec_cn(cn,17))/(1-cn_omega);
            iec_cn(cn,11)=log10(1-cn_omega)-enerate/iec_cn(cn,14);
            iec_cn(cn,24)=iec_cn(cn,11);
            iec_cn(cn,25)=0;
            iec_cn(cn,26)=log10(1-cn_omega);
            iec_cn(cn,27)=enerate;
        end
        
        % setting up TFP in DICE
        if fudanccm_exo==3
            tfprate = 0.016*exp(-0.001*(realtime(t,1)-2019));
            if dicetfp>=1
                if t>121 && cn>1
                    tfprate = (1-0.3)* cn_growth *exp(-0.001*(realtime(t,1)-2019));
                    % Burke 2015 function for TFP
                    if iec_cn(cn,29)==1
                        temp_tfp = burkecoef(1)*(cntemp^2-cntemp0^2)-tempopt*burkecoef(1)*2*(cntemp-cntemp0);
                    else
                        temp_tfp = burkecoef(2)*(cntemp^2-cntemp0^2)-tempopt*burkecoef(1)*2*(cntemp-cntemp0);
                    end
                    tfprate = tfprate - max(-tfprate, temp_tfp); % avoiding bad value in Monte Carlo
                    if dicetfp==1
                        tfprate = max(-0.02, tfprate);
                    else
                        tfprate = max(0, tfprate);
                    end
                end
            end
            cn_econ0(3) = cn_econ0(3) * (1+tfprate)^realtime(t,2); % TFP in DICE
        end
        
        %Economy
        S(t+1,1:23,cn+1)=econdyn(t, cn_labor, cn_econ0, cn_fracinv, iec_cn(cn,:), cn_omega, cn_invest, LR, D, inertia, deff, switcher);
        
        %Considering learning in DICE for dicetfp=3 only
        if fudanccm_exo==3 && dicetfp<3
            if t<=49
                S(t,4,cn+1) = 550; % The initial price  550 USD (t CO2)-1, which decreases at a constant rate of 0.5% y-1
            end
            if dicetfp<3
                S(t+1,4,cn+1) = S(t,4,cn+1) * (1-0.005)^realtime(t,2);
                S(t+1,19,cn+1) = S(t+1,4,cn+1) * S(t+1,18,cn+1)^(theta2-1);
            end
        end
        
        %Population
        S(t+1,35,cn+1)=S(t,35,cn+1)*L(t+1)/L(t); % labor
        
        %Utility
        if (t+1)>=121 && (t+1)<=400
            S(t+1,36,cn+1) = (((1-rsav)* S(t+1,7,cn+1)/S(t+1,35,cn+1)*1000)^(1-elasmu)/(1-elasmu)) * S(t+1,35,cn+1) /1000  * (1-rour*0.001)^(t-120); % utility discounted to 2025  
        end
    end
    
    %Global/regional technology diffusion
    if fudanccm_exo==1 && techdiff>0 && realtime(t,1)>2020
        ZCE = zeros(5,2);
        for cn=2:cn_num
            zz=cou_iform(cndata(cn,1),4);
            ZCE(zz,1) = ZCE(zz,1) + S(t,21,cn+1); % green energy PJ
            ZCE(zz,2) = ZCE(zz,2) + S(t+1,21,cn+1); % green energy PJ
            ZCE(5,1) = ZCE(5,1) + S(t,21,cn+1); % green energy PJ
            ZCE(5,2) = ZCE(5,2) + S(t+1,21,cn+1); % green energy PJ
        end
        for zz=1:5
            ZCE(zz,1) = max(72000,ZCE(zz,1)); % energy PJ
            ZCE(zz,2) = max(ZCE(zz,1),ZCE(zz,2)); % energy PJ
        end
        if ZCE0(5,1)==0
            ZCE0=ZCE; % initializing
        end
        for cn=2:cn_num
            S(t+1,4,cn+1) = S(t,4,cn+1) * (min(1+0.05*realtime(t,2), max(1e4, S(t+1,12,cn+1)) / max(1e4, S(t,12,cn+1))))^(theta2-1); % Effect of energy expansion
            if techdiff==1
                % 3 global technology diffusion
                S(t+1,4,cn+1) = S(t+1,4,cn+1)*(max(1+0.01*realtime(t,2), ZCE(5,2)/ZCE(5,1) ))^(log2(1-min(LR*2,max(LR,LR+log ( ZCE(5,1) /3600/50) / log(4) * 0.2))));
            elseif techdiff==2
                % 4 regional technology diffusion
                ZCEnergy0 = max(72000,ZCE(zz,1)-ZCE0(zz,1)+ZCE0(5,1));
                ZCEnergy1 = max(ZCEnergy0,ZCE(zz,2)-ZCE0(zz,1)+ZCE0(5,1));
                S(t+1,4,cn+1) = S(t+1,4,cn+1)*(max(1+0.01*realtime(t,2), ZCEnergy1/ZCEnergy0 ))^(log2(1-min(LR*2,max(LR,LR+log ( ZCEnergy0 /3600/50) / log(4) * 0.2))));
            end
            S(t+1,19,cn+1) = S(t+1,4,cn+1) * S(t+1,18,cn+1)^(theta2-1); % carbon price
        end
    end
    if fudanccm_exo==1 && techdiff>0 && realtime(t,1)<=2020
        for cn=2:cn_num
            S(t+1,4,cn+1) = 550;
            S(t+1,19,cn+1) = S(t+1,4,cn+1) * S(t+1,18,cn+1)^(theta2-1); % carbon price
        end
    end
    
    %Nuclear Energy
    totalemico2=sum(S(t+1,20,3:(cn_num+1)),3); % total emissions
    if nuclear_tip>0 && totalemico2>0 && t>=121 && realtime(t+1,1)>abtcn(cn,7)
        % initializing data before 2025
        if t==121
            for t_nu=1:121
                S(t_nu,38,3:(cn_num+1)) = cnenergy(2:cn_num,min(t_nu,49),6) * 0.0928; % quad Btu = 0.29 PWh * heat rate of 32.63% = 0.0928
            end
        end
        
        % total gdp
        totalgdpt=sum(S(t,13,3:(cn_num+1)),3);
        
        % cumulative nuclear power by region
        nuclearregion = zeros(5,1);
        for cn=2:cn_num
            zz=cou_iform(cndata(cn,1),4);
            nuclearregion(zz,1) = nuclearregion(zz,1) + sum(S(121:t,38,cn+1))+ (sum(cnenergy(cn,1:49,6),2)+cnenergy(cn,49,6)*5) * 0.0928; % cumulative nuclear energy PWh
        end
        nuclearregion(5,1) = sum( nuclearregion(1:4,1),1);
        
        for cn=2:cn_num
            if nuclear_tech==6
                nuclear_LF=1-1/20;
            else
                nuclear_LF=1-1/60;
            end
            % CO2 emissions abated by nuclear power with a lifetime of 60 years 1-1/60=0.9833
            emifrac=S(t+1,20,cn+1)/totalemico2;  % fraction of CO2 emissions in global total
            if (realtime(t+1,1)-abtcn(cn,7))<=abtcn(cn,9)
                newnuclear = abtcn(cn,8) * emifrac * min(1,max(0,(realtime(t+1,1)-abtcn(cn,7))/abtcn(cn,9))); % new nuclear power TW deployed in a country
                newnuclear2 = abtcn(cn,8) * emifrac / abtcn(cn,9); % add nuclear to reach the TW
            elseif (realtime(t+1,1)-abtcn(cn,7))>abtcn(cn,9)
                newnuclear = abtcn(cn,8) * emifrac * min(1,max(0,(realtime(t+1,1)-abtcn(cn,7))/abtcn(cn,9))); % new nuclear power TW deployed in a country
                if S(t+1,20,cn+1)>=(newnuclear * 8.76 * 0.8433 * 0.85)
                    newnuclear2 = newnuclear*(1-nuclear_LF); % add nuclear to sustain the TW
                else
                    newnuclear = newnuclear*nuclear_LF;  % let the nuclear TW to decline
                    newnuclear2 = 0;
                end
            end
            S(t+1,37,cn+1) = S(121,38,cn+1) / 6.132 * 0.9833^(realtime(t+1,1)-abtcn(cn,7)) + newnuclear; % Total nuclear power TW based on the capacity factor of 70% for early years(Koomey,2007) 365*24*0.7/1000=6.132 PWh/TW
            S(t+1,40,cn+1) = min(0.85, S(t+1,20,cn+1) / max(0.000001, newnuclear * 8.76 * 0.8433)); % Real capacity factor based on 365*24/1000=8.76 and FF CO2 emission factor of 0.8433 t CO2 / MWh
            S(t+1,38,cn+1) = S(121,38,cn+1) * 0.9833^(realtime(t+1,1)-2025) + newnuclear * 8.76 * S(t+1,40,cn+1); % Nuclear energy PWh based on a capacity factor of 85% for the future (Koomey,2007)
            S(t+1,39,cn+1) = newnuclear * 8.76 * 0.8433 * S(t+1,40,cn+1); % Abated CO2 emissions by new nuclear GtCO2/yr based on FF CO2 emission factor of 0.8433 t CO2 / MWh
            
            % accumulation of nuclear power PWh after 1 PWh minimum for global technological baseline
            % nuclear_tech: 0 no tech advances; 1 national tech advances (default); 2 global cooperation; 3 OECD/REF/ALM/ASIA cooperation; 4 no accident
            if nuclear_tech==0
                cum_nuclear = max(1,(sum(cnenergy(cn,1:49,6),2)+cnenergy(cn,49,6)*5) * 0.0928); % nuclear PWh from 1980 to 2025
            elseif nuclear_tech==2
                cum_nuclear = nuclearregion(5,1);
            elseif nuclear_tech==3
                cum_nuclear = nuclearregion(cou_iform(cndata(cn,1),4),1);
            else
                cum_nuclear = sum(S(121:t,38,cn+1)) + max(1,(sum(cnenergy(cn,1:49,6),2)+cnenergy(cn,49,6)*5) * 0.0928); % national tech advances (default)
            end
            
            % cumulative nuclear power PWh that will have a new accident
            log10_nuclear = 0.707 * log10(cum_nuclear) - 1.585;             
            if nuclear_tip==2
                log10_nuclear = log10_nuclear + (-0.092 * log10(cum_nuclear) + 0.732) * randnumc(1,6); % nuclear expectable of a new accident: normal distribution
%                 idddd = idddd + 1; idddd = mod(idddd,1000);
%                 log10_nuclear = log10_nuclear + (-0.125 * log10(cum_nuclear) + 0.779) * randnumc(1,end-idddd); % nuclear expectable of a new accident: normal distribution
            end
            cri_nuclear=10^log10_nuclear;
            id_nuclear(cn-1,2)=id_nuclear(cn-1,2)+S(t+1,38,cn+1); % PWh cumlative nuclear since last accident
            
            % nuclear accident
            Dnuclear=0;
            if nuclear_tech~=4 && id_nuclear(cn-1,2)>cri_nuclear
                S(t+1,41,cn+1) = floor( id_nuclear(cn-1,2)/cri_nuclear); % number of accident
                id_nuclear(cn-1,2) = id_nuclear(cn-1,2) - S(t+1,41,cn+1) * cri_nuclear;
                for nui=1:S(t+1,41,cn+1)
                    if nuclear_tip==1
                        nulevel=5;
                    elseif nuclear_tip==2
                        % id of nuclear accident in this country
                        if id_nuclear(cn-1,1)==0
                            id_nuclear(cn-1,1) = cn*3-5+6;
                        elseif id_nuclear(cn-1,1)==3156
                            id_nuclear(cn-1,1)=7;
                        else
                            id_nuclear(cn-1,1)=id_nuclear(cn-1,1)+1;
                        end
                        nulevel = randnumc(1,id_nuclear(cn-1,1)); % level of nuclear accident 1-7
                    end
                    % damage of nuclear accident % to GDP
                    if nuclear_acc==1
                        Dnuclear=Dnuclear+10^(nulevel*0.684-6.268)*S(t,13,cn+1); % trillion $
                    elseif nuclear_acc==2
                        Dnuclear=Dnuclear+10^(nulevel*0.901-6.836)*S(t,13,cn+1); % trillion $
                    elseif nuclear_acc==3
                        Dnuclear=Dnuclear+10^(nulevel*0.534-5.619); % trillion $
                    elseif nuclear_acc==4
                        if nulevel<5
                        Dnuclear=Dnuclear+10^(nulevel*0.438-2.329-3);  
                        else
                        Dnuclear=Dnuclear+10^(nulevel*1.181-6.130-3); 
                        end% trillion $
                        
                    end
                    if nulevel>4
                        S(t+1,44,cn+1)=S(t+1,44,cn+1)+1; % number of accident level>=5
                    end
                end
            end
            
            % nuclear accident damage as a percentage to GDP
            if Dnuclear<=(S(t,13,cn+1)*0.5)
                S(t+1,42,cn+1) = min(0.9, S(t+1,42,cn+1)+Dnuclear/S(t,13,cn+1)); % as a percentage to GDP for this country only
            else
                globaldamage = (Dnuclear-S(t,13,cn+1)*0.5) / totalgdpt;
                S(t+1,42,3:(cn_num+1)) = min(0.9, S(t+1,42,3:(cn_num+1)) + globaldamage); % as a percentage to GDP for all countries                    
                S(t+1,42,cn+1) = min(0.9, S(t+1,42,cn+1) + 0.5 - globaldamage);  % 50% of GDP lost for this country only
            end
            
            % cost of deploying nuclear power trillion $
            if nuclear_tech==5
                S(t+1,43,cn+1) = newnuclear2 * 2 + S(t+1,38,cn+1) * 0.0158; % construction costs of 2000 $/kW and 0.0055 $/kWh for uranium cost (51% in operational costs)
            else
                S(t+1,43,cn+1) = newnuclear2 * 4 + S(t+1,38,cn+1) * 0.0158; % default: construction costs of 4000 $/kW and 0.0055 $/kWh for uranium cost (51% in operational costs)
            end
            
            % damage to the economy
            S(t+1,5,cn+1)=S(t+1,5,cn+1)+S(t+1,43,cn+1)/S(t,13,cn+1); % nuclear accelerate EUE            
%             S(t+1,6,cn+1)=S(t+1,6,cn+1)+S(t+1,42,cn+1); % nuclear decelerate ENE
%             S(t+1,7,cn+1)=max(S(t+1,7,cn+1)*0.1, S(t+1,7,cn+1)-S(t,13,cn+1)*S(t+1,42,cn+1)-S(t+1,43,cn+1)); % reduce net output for both consumption and investment
            S(t+1,6,3:(cn_num+1))=S(t+1,6,3:(cn_num+1)) + Dnuclear / totalgdpt; % nuclear decelerate ENE            
            S(t+1,13,cn+1)=S(t+1,13,cn+1)-S(t+1,43,cn+1); % GDP reduced by the costs of deploying nuclear power
            S(t+1,13,3:(cn_num+1))=S(t+1,13,3:(cn_num+1)) .* (1-Dnuclear / totalgdpt); % GDP in all countries reduced by the damage of nuclear accidents
        end
    end
    
    %Considering fast energy production growth in developing countries
    if energy_growth>0 && t>=121
        for cn=2:cn_num
            pce2025 = max(1, energy_growth / (S(121,12,cn+1) / 3600 * 0.3263 * 1000 / S(121,35,cn+1))); % per capita energy in MWh/cap
            if pce2025>1
                S(t+1,2,cn+1) = S(t,2,cn+1) * (S(121,2,cn+1) * pce2025 / S(t,2,cn+1))^(0.01 * realtime(i,2)) ;
            end
        end
    end
    
    %Using the MACC by project
    if mac_curves==1 && t==120
        amac=zeros(cn_num*2-2,5);
        % assume a baseline MACC (do not change this part)
        for cn=1:(cn_num-1)
            % a set of expensive abatement
            amac(cn,1)=300; % lowest mac $/tCO2
            amac(cn,2)=0; % lowest emission abatement Gt CO2
            amac(cn,3)=300*3; % highest mac $/tCO2
            amac(cn,4)=S(121,20,cn+2)*2; % highest emission abatement Gt CO2
            amac(cn,5)=cn; % id of country
        end
        % create an assumed MACC (this part is to be updated)
        for cn=1:(cn_num-1)
            amac(cn+cn_num-1,1)=0; % lowest mac $/tCO2
            amac(cn+cn_num-1,2)=0; % lowest emission abatement Gt CO2
            amac(cn+cn_num-1,3)=300; % highest mac $/tCO2
            amac(cn+cn_num-1,4)=S(121,20,cn+2); % highest emission abatement Gt CO2
            amac(cn+cn_num-1,5)=cn; % id of country
        end
        % MACC for 134 countries
        globalmac=globalmac_curve(amac, cn_num-1);
    end
    
    %Climate Investment Transfer
    if tempw>1 && invest_trans(3)>0
        iec_cn(:,17)=0; % EUE rate calibrated to 0 for no mitigation
        if t<120
            % initializing
            for cn=1:(cn_num-1)
                abtcn(cn+1,4)=exp(log(100/ (S(t+1,7,cn+2) * 1000 / S(t+1,35,cn+2))) /100) - 1; % per capita gdp reaches 100k$ after 100 years
            end
        else
            % determine the ratio of investment after exceeding a threshold to meet the target of global climate investment
            cn_investt0=0; % global climate investment using 100% GDP
            for cn=1:(cn_num-1)
                cn_pcg = max(0, S(t+1,7,cn+2) * 1000 / S(t+1,35,cn+2) - invest_trans(2)); % per capita gdp k$
                cn_investt0 = cn_investt0 + cn_pcg * S(t+1,35,cn+2) /1000; % climate investment trillion $
            end
            gabtcost = sum(globalmac(:,1:(cn_num-1),2),2); % global abatement cost as a function of mac trillion $
            invest_ratio = min(0.5, sum(S(t+1,7,3:(cn_num+1)),3) * invest_trans(3) / cn_investt0); % fraction of GDP used as climate investment
            for cn=1:(cn_num-1)
                cn_pci = max(0, S(t+1,7,cn+2) * 1000 / S(t+1,35,cn+2) - invest_trans(2))  * invest_ratio; % per capita climate investment k$
                S(t+1,45,cn+2) = cn_pci * S(t+1,35,cn+2) /1000; % climate investment trillion $
            end
            cn_investt1 = sum(S(t+1,45,3:(cn_num+1)),3); % global climate investment trillion $
            idx_invest=find(gabtcost<=cn_investt1);
            gcp=min(size(gabtcost,1)-1, idx_invest(end)); % global carbon price for all countries based on 100% of climate investment transfer
            if (gabtcost(gcp+1) - gabtcost(gcp))>0
                fextrainvest = (cn_investt1 - gabtcost(gcp))/(gabtcost(gcp+1) - gabtcost(gcp)); % extra emission abatement fraction
            else
                fextrainvest = 0;
            end
            
            % matrix of climate investment transfer between countries
            for cn=1:(cn_num-1)
                cn_investt2 = globalmac(gcp,cn,2)+(globalmac(gcp+1,cn,2)-globalmac(gcp,cn,2))*fextrainvest; % optimal climate investment trillion $
                S(t+1,46,cn+2) = S(t+1,45,cn+2) + (cn_investt2-S(t+1,45,cn+2)) * invest_trans(1); % climate investment used in a country trillion $
            end
            
            % climate investment transfer from region A to region B
            if invest_trans(1)>0
                % source of climate investment transfer
                invest_flux = zeros(cn_num,2);
                for cn=1:(cn_num-1)
                    if S(t+1,45,cn+2)>S(t+1,46,cn+2)
                        invest_flux(cn,1) = S(t+1,45,cn+2) - S(t+1,46,cn+2); % climate investment transfer from a country
                    else
                        invest_flux(cn,2) = S(t+1,46,cn+2) - S(t+1,45,cn+2); % climate investment transfer to a country
                    end
                end
                invest_flux(cn_num,1:2)=sum(invest_flux(1:(cn_num-1),1:2),1);
                % climate investment transfer from cn2 to cn1
                invest_matrix = zeros(cn_num-1,cn_num-1); % positive for getting money
                for cn1=1:(cn_num-1)
                    if invest_flux(cn1,2)>0
                        for cn2=1:(cn_num-1)
                            if invest_trans(4)==1
                                % block climate investment transfer to a country / region
                                if cou_iform(cndata(cn1+1,1),3)==invest_trans(6) && cou_iform(cndata(cn2+1,1),3)==invest_trans(5)
                                    continue;
                                end
                            elseif invest_trans(4)==2
                                % open climate investment transfer to a country / region
                                if cou_iform(cndata(cn1+1,1),3)~=invest_trans(6) || cou_iform(cndata(cn2+1,1),3)~=invest_trans(5)
                                    continue;
                                end
                            end
                            invest_matrix(cn1,cn2)=invest_flux(cn1,2)*invest_flux(cn2,1)/invest_flux(cn_num,1); % transfer from cn2 to cn1
                            invest_matrix(cn2,cn1)=-invest_matrix(cn1,cn2); % transfer from cn1 to cn2
                        end
                    end
                end
                for cn1=1:(cn_num-1)
                    S(t+1,46,cn1+2) = S(t+1,45,cn1+2) + sum(invest_matrix(cn1,:),2);
                    for cn2=1:(cn_num-1)
                        zz=cou_iform(cndata(cn2+1,1),3);
                        S(t+1,48+zz,cn1+2) = S(t+1,48+zz,cn1+2) + invest_matrix(cn1,cn2);
                    end
                end
            end
            
            % determine the national emission abatement
            for cn=1:(cn_num-1)
                if S(t+1,46,cn+2)>0
                    idx_investcn=find(globalmac(:,cn,2)<=S(t+1,46,cn+2)); % seek for national carbon price
                    if size(idx_investcn,1)>1
                        if idx_investcn(end)==size(globalmac,1)
                            % all emission abatement projects are realized
                            S(t+1,46,cn+2) = globalmac(end,cn,2); % climate investment used in a country trillion $
                            S(t+1,47,cn+2) = globalmac(end,cn,1); % emission abatement in a country Gt CO2
                            S(t+1,48,cn+2) = size(globalmac,1)+1; % national carbon price $/tCO2
                        else
                            gcpcn=idx_investcn(end); % national carbon price $/tCO2
                            if (globalmac(gcpcn+1,cn,2) - globalmac(gcpcn,cn,2))>0
                                fextrainvestcn = (S(t+1,46,cn+2) - globalmac(gcpcn,cn,2)) / (globalmac(gcpcn+1,cn,2) - globalmac(gcpcn,cn,2)); % extra emission abatement fraction
                            else
                                fextrainvestcn = 0;
                            end
                            S(t+1,47,cn+2) = globalmac(gcpcn,cn,1) + (globalmac(gcpcn+1,cn,1) - globalmac(gcpcn,cn,1)) * fextrainvestcn; % emission abatement in a country Gt CO2
                            S(t+1,48,cn+2) = gcpcn+fextrainvestcn-1; %  national carbon price $/tCO2
                        end
                    end
                end
            end
            investmentchange=sum(S(t+1,46,3:(cn_num+1)),3)/sum(S(t+1,45,3:(cn_num+1)),3);
            S(t+1,45,3:(cn_num+1))=S(t+1,45,3:(cn_num+1)).*investmentchange;
            
            % update global mac curve for the next year
            for cn=1:(cn_num-1)
                if S(t+1,47,cn+2)>(S(t+1,20,cn+2)*0.1)
                    idx_learning=find(amac(cn_num:end,5)==cn);
                    if size(idx_learning,1)>0
                        amac(cn_num-1+idx_learning,1:2)=amac(cn_num-1+idx_learning,1:2).*(1-0.005)^realtime(t,2); % MAC reduced by mitigation
                    end
                end
            end
            globalmac = globalmac_curve(amac, cn_num-1);
            
            % damage to the economy
            for cn=1:(cn_num-1)
                S(t+1,5,cn+2) = S(t+1,5,cn+2) + S(t+1,46,cn+2)/2 /S(t+1,7,cn+2);
                S(t+1,6,cn+2) = S(t+1,6,cn+2) + ( S(t+1,45,cn+2) - S(t+1,46,cn+2)/2 )/S(t+1,7,cn+2);
                S(t+1,13,cn+2) = S(t+1,13,cn+2) - S(t+1,45,cn+2);
                if S(t+1,3,cn+2) < S(t,3,cn+2)
                    S(t+1,3,cn+2) = S(t,3,cn+2); % ensuring ENE growth
                end
            end
        end
    end
    
    %Negative emissions
    if realtime(t,1)<neg_year
        emi_neg0 = emi_neg * max(0,(realtime(t,1)-2025)/(neg_year-2025)); % deploying negative emissions since 2025
    else
        
        emi_neg0 = emi_neg; % deploying emi_neg after 2025
    end

    if tempw>1 && (emi_neg0+max(abtcn(2:cn_num,5),[],1)+max(abtcn(2:cn_num,6),[],1))>0
        emico2_pre=0; % CO2 emission Gt CO2 without abatement
        for cn=2:cn_num
            if fudanccm_exo==3
                emico2_pre=emico2_pre+S(t+1,13,cn+1) * S(t+1,1,cn+1); % DICE
            else
                emico2_pre=emico2_pre+S(t+1,12,cn+1) * EFco2(min(45,t)) * (1-Egreen(min(45,t+1),8)); % FUDAM
            end
        end
        fract_neg=emi_neg0/emico2_pre; % global fraction of negative emissions relative to the fossil fuel emissions without abatement
        % for tipping with nuclear power
        if nuclear_tip>0 && realtime(t,1)>=2200
            fract_neg=fract_neg*2;
        end
        
        for cn=2:cn_num
            if fudanccm_exo==3
                fract_cnneg = min(3,abtcn(cn,6) / (S(t+1,13,cn+1) * S(t+1,1,cn+1))) + abtcn(cn,5);
                S(t+1,20,cn+1) = S(t+1,20,cn+1) - S(t+1,13,cn+1) * S(t+1,1,cn+1) * S(t+1,18,cn+1) * (fract_neg + fract_cnneg);
            else
                fract_cnneg = min(3,abtcn(cn,6) / (S(t+1,12,cn+1) * EFco2(min(45,t)) * (1-Egreen(min(45,t+1),8)))) + abtcn(cn,5);
                S(t+1,20,cn+1) = S(t+1,20,cn+1) - S(t+1,12,cn+1) * EFco2(min(45,t)) * (1-Egreen(min(45,t+1),8)) * S(t+1,18,cn+1) * (fract_neg + fract_cnneg);
            end
            S(t+1,18,cn+1) = S(t+1,18,cn+1)*(1+fract_neg+fract_cnneg); % Fraction of emission abatement
            S(t+1,19,cn+1) = S(t+1,4,cn+1)*S(t+1,18,cn+1)^(theta2-1); % Carbon price
        end
    end

    %CO2 emission Gt CO2
    emico2=sum(S(t+1,20,3:(cn_num+1)),3);
    
    %Additional emissions after a tipping function according to Cai 2016 NCC
    if tipping>0
        for tip_i=2:5
            if tip_prob(tip_i,2)==0
                if (t-tip_prob(tip_i,4))<=tip_para(tip_i,2) || tip_i==5
                    emico2 = emico2 + tip_para(tip_i,4) * 3.666; % Gt CO2
                end
            end
        end
    end
    
    % Estimate scc by adding emissions
    if floor(realtime(t,1))==scc(1) && scc(2)~=0
        emico2=emico2+scc(2); % Gt CO2
    end
    
    % CO2 emission abated by nuclear power
    if nuclear_tip>0 && t>=121
        emico2=emico2-sum(S(t+1,39,3:(cn_num+1)),3);
    end
    
    % CO2 emission abated by international climate finance
    if invest_trans(3)>0 && tempw>1
        emico2=emico2-sum(S(t+1,47,3:(cn_num+1)),3);
    end
    
    % Climate module
    clim2=climdyn(t, clim1, FFlux, emico2);
    t=t+1;
    econ1=econ2;
    clim1=clim2;
    S(t,1:23,1)=econ2(1,1:23);
    S(t,24:34,1)=clim2(1,1:11);
    S(t,35,1)=L(t); % labor
    
end

for cn=1:(cn_num+1)
    S(:,1,cn) = S(:,1,cn) * 3600; % EUE $/KJ -> $/kWh
    S(:,2,cn) = S(:,2,cn) / 3600; % EPE PJ/(t$)^0.3 / (billion cap)^0.7 -> PWh/(t$)^0.3 / (billion cap)^0.7
    S(:,12,cn) = S(:,12,cn) / 3600 * 0.3263; % energy PJ -> PWh based on a heat rate of 32.63%
    S(:,21,cn) = S(:,21,cn) / 3600 * 0.3263; % cumulative green energy PJ -> PWh based on a heat rate of 32.63%
end

tend=endyear-2025+120;
npv=sum(sum(S(121:tend,36,3:(cn_num+1)),3),1); % global total NPV of utility

end
