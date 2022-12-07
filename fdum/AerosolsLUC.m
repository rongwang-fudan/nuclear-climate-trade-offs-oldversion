% Author: Rong Wang // contact rongwang@fudan.edu.cn //
% Date: 2021.2.10

global  carbonbudget20 Fex Eland realtime

% carbonbudget20
% 1 year; 2 fossil emi GtC; 3 LUC emi GtC; 4 Atm growth GtC; 5 Ocean sink GtC; 6 Land sink GtC; 7 Cement carbonation
% Budget = 2 + 3 - 4 - 5 - 6 - 7

% Time horizon
T = size(realtime,1);

% Radiative forcing by 1 CH4, 2 N2O, 3 CFCs, 4 aerosol
Fex = zeros(T,4);
% CH4
% 0.1 in 1900, to 0.25 in 1950, and to 0.4 in 1990, to 0.5 in 2100, and to 0.1 in 2500
% Figure 8.6 from IPCC-AR5-Chapter 8
% lifetime_CH4=12; % http://chartsbin.com/view/2407
for i=1:T
    if realtime(i,1)<=1900
        Fex(i,1) = 0.1;
    elseif realtime(i,1)<=1950
        Fex(i,1) = 0.1+(realtime(i,1)-1900)*0.15/50;
    elseif realtime(i,1)<=1990
        Fex(i,1) = 0.25+(realtime(i,1)-1950)*0.15/40;
    elseif realtime(i,1)<=2100
        Fex(i,1) = 0.4+(realtime(i,1)-1990)*0.1/110;
    else
        Fex(i,1) = 0.5-(realtime(i,1)-2100)*0.4/400;
    end
end

% N2O
% 0.01 in 1900, to 0.06 in 1950, to 0.18 in 2010, to 0.2 in 2100, and to 0 in 2050
% Figure 8.6 from IPCC-AR5-Chapter 8
% halflife_N2O=114 * 0.693; % http://chartsbin.com/view/2407
for i=1:T
    if realtime(i,1)<=1900
        Fex(i,2) = 0.01;
    elseif realtime(i,1)<=1950
        Fex(i,2) = 0.01+(realtime(i,1)-1900)*0.05/50;
    elseif realtime(i,1)<=2010
        Fex(i,2) = 0.06+(realtime(i,1)-1950)*0.12/60;
    elseif realtime(i,1)<=2100
        Fex(i,2) = 0.18+(realtime(i,1)-2010)*0.02/90;
    else
        Fex(i,2) = 0.2-(realtime(i,1)-2100)*0.19/400;
    end
end

% CFCs
% 0 in 1900, to 0.01 in 1960, to 0.3 in 1990, to 0.4 in 2100, and to 0 in 2050
% Figure 8.6 from IPCC-AR5-Chapter 8
% halflife_CFCs=25 * 0.693; % http://chartsbin.com/view/2407
for i=1:T
    if realtime(i,1)<=1900
        Fex(i,3) = 0;
    elseif realtime(i,1)<=1960
        Fex(i,3) = (realtime(i,1)-1900)*0.01/60;
    elseif realtime(i,1)<=1990
        Fex(i,3) = 0.01+(realtime(i,1)-1960)*0.29/30;
    elseif realtime(i,1)<=2100
        Fex(i,3) = 0.3+(realtime(i,1)-1990)*0.1/110;
    else
        Fex(i,3) = 0.4-(realtime(i,1)-2100)*0.4/400;
    end
end

% AEROSOL
% -0.08 in 1900, to -0.1 in 1950, to -0.3 in 1970, and to -0.4 in 1990, to -0.5 in 2010, to -0.1 in 2100
% Figure 8.8 from IPCC-AR5-Chapter 8
for i=1:T
    if realtime(i,1)<=1900
        Fex(i,4) = -0.08;
    elseif realtime(i,1)<=1950
        Fex(i,4) = -0.08-(realtime(i,1)-1900)*0.02/50;
    elseif realtime(i,1)<=1970
        Fex(i,4) = -0.1-(realtime(i,1)-1950)*0.2/20;
    elseif realtime(i,1)<=1990
        Fex(i,4) = -0.3-(realtime(i,1)-1970)*0.1/20;
    elseif realtime(i,1)<=2010
        Fex(i,4) = -0.4-(realtime(i,1)-1990)*0.1/20;
    elseif realtime(i,1)<=2100
        Fex(i,4) = -0.5+(realtime(i,1)-2010)*0.4/90;
    else
        Fex(i,4) = -0.1;
    end
end

% Carbon emissions from land use change
% Carbonbudget20
%   1 year; 2 fossil emi GtC; 3 LUC emi GtC; 4 Atm growth GtC; 5 Ocean sink GtC; 6 Land sink GtC; 7 Cement carbonation
Eland = zeros(T,1);
for i=1:T
    if i<=size(carbonbudget20,1)
        Eland(i,1) = carbonbudget20(i,3);
    else
        Eland(i,1) = 0.2 + (Eland(i-1,1)-0.2) * (1-0.023)^realtime(i,2); % decline at 2.3% yr-1
    end
end
