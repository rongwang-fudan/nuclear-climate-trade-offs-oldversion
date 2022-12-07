% Author: Rong Wang // contact rongwang@fudan.edu.cn //
% Date: 2021.7.28
% Calibtration of efficiency change during COVID-19
%    delta_eff: changes in 1 EUE, 2 EPE, and 3 ENE
% Output (Y) is sensitivty to ENE; Energy (E) is sensitivty to EUE; and Price (pe) and energy (E) are sensitivty to EPE

function [ output_covid, delta_eff] = Calibration_COVID( FFlux, L, iec, calrsav, LR, switcher, covidyear, displays, sensitivity )

global S0

if sensitivity==1
    iec(1,2) = iec(1,2) * (1 + iec(1,3)*1.96);
    iec(1,4) = iec(1,4) * (1 + iec(1,5)*1.96);
    iec(1,12) = iec(1,12) * (1 + iec(1,13)*1.96);
    iec(1,14) = iec(1,14) * (1 + iec(1,15)*1.96);
elseif sensitivity==2
    iec(1,2) = iec(1,2) * (1 - iec(1,3)*1.96);
    iec(1,4) = iec(1,4) * (1 - iec(1,5)*1.96);
    iec(1,12) = iec(1,12) * (1 - iec(1,13)*1.96);
    iec(1,14) = iec(1,14) * (1 - iec(1,15)*1.96);
end

% Finding the optimal delta_eff to achieve the best agreement with observations
deff = [0, 0, 0.05, 0.05, 0.08, 0, 0, 0, 0, -0.01, -0.02;
    0, 0, 0, 0, 0.01, 0.02, 0.01, 0.01, 0.01, 0, 0;
    -0.06, -0.06, -0.02, 0.02, 0.04, 0, 0, 0, 0, 0, 0];

[S2, LFmin] = covid( FFlux, L, iec, calrsav, LR, switcher, covidyear, deff );
deff2 = deff;

%Calibration of Y by ENE
deffnew=deff;
deffnew(3,:)=deffnew(3,:)*1.1;
[ss, LF] = covid( FFlux, L, iec, calrsav, LR, switcher, covidyear, deffnew );
if LF(end,2)<=LFmin(end,2)
    while LF(end,2)<=LFmin(end,2) && (deff(3,5)/deff2(3,5))<5
        if displays==1
            display('ENE+1');
        end
        S2=ss;
        deff=deffnew;
        LFmin=LF;
        % try to use new parameters
        deffnew(3,:)=deffnew(3,:)*1.1;
        [ss, LF] = covid( FFlux, L, iec, calrsav, LR, switcher, covidyear, deffnew );
    end
else
    deffnew=deff;
    deffnew(3,:)=deffnew(3,:)*0.9;
    [ss, LF] = covid( FFlux, L, iec, calrsav, LR, switcher, covidyear, deffnew );
    while LF(end,2)<=LFmin(end,2) && (deff(3,5)/deff2(3,5))>0.2
        if displays==1
            display('ENE-1');
        end
        S2=ss;
        deff=deffnew;
        LFmin=LF;
        % try to use new parameters
        deffnew(3,:)=deffnew(3,:)*0.9;
        [ss, LF] = covid( FFlux, L, iec, calrsav, LR, switcher, covidyear, deffnew );
    end
end

%Calibration of E by EUE
deffnew=deff;
deffnew(1,:)=deffnew(1,:)*1.1;
[ss, LF] = covid( FFlux, L, iec, calrsav, LR, switcher, covidyear, deffnew );
if LF(end,1)<=LFmin(end,1)
    while LF(end,1)<=LFmin(end,1) && (deff(1,5)/deff2(1,5))<5
        if displays==1
            display('EUE+1');
        end
        S2=ss;
        deff=deffnew;
        LFmin=LF;
        % try to use new parameters
        deffnew(1,:)=deffnew(1,:)*1.1;
        [ss, LF] = covid( FFlux, L, iec, calrsav, LR, switcher, covidyear, deffnew );
    end
else
    deffnew=deff;
    deffnew(1,:)=deffnew(1,:)*0.9;
    [ss, LF] = covid( FFlux, L, iec, calrsav, LR, switcher, covidyear, deffnew );
    while LF(end,1)<=LFmin(end,1) && (deff(1,5)/deff2(1,5))>0.2
        if displays==1
            display('EUE-1');
        end
        S2=ss;
        deff=deffnew;
        LFmin=LF;
        % try to use new parameters
        deffnew(1,:)=deffnew(1,:)*0.9;
        [ss, LF] = covid( FFlux, L, iec, calrsav, LR, switcher, covidyear, deffnew );
    end
end

%Calibration of pe by EPE
deffnew=deff;
deffnew(2,:)=deffnew(2,:)*1.1;
[ss, LF] = covid( FFlux, L, iec, calrsav, LR, switcher, covidyear, deffnew );
if LF(end,3)<=LFmin(end,3) && LF(end,1)<=LFmin(end,1)
    while LF(end,3)<=LFmin(end,3) && LF(end,1)<=LFmin(end,1) && (deff(2,5)/deff2(2,5))<5
        if displays==1
            display('EPE+1');
        end
        S2=ss;
        deff=deffnew;
        LFmin=LF;
        % try to use new parameters
        deffnew(2,:)=deffnew(2,:)*1.1;
        [ss, LF] = covid( FFlux, L, iec, calrsav, LR, switcher, covidyear, deffnew );
    end
else
    deffnew=deff;
    deffnew(2,:)=deffnew(2,:)*0.9;
    [ss, LF] = covid( FFlux, L, iec, calrsav, LR, switcher, covidyear, deffnew );
    while LF(end,3)<=LFmin(end,3) && LF(end,1)<=LFmin(end,1) && (deff(2,5)/deff2(2,5))>0.2
        if displays==1
            display('EPE-1');
        end
        S2=ss;
        deff=deffnew;
        LFmin=LF;
        % try to use new parameters
        deffnew(2,:)=deffnew(2,:)*0.9;
        [ss, LF] = covid( FFlux, L, iec, calrsav, LR, switcher, covidyear, deffnew );
    end
end

delta_eff=deff;
S0(45:126,1:34) = S2(1:82,1:34);

output_covid = S2(:,1:24);
output_covid(:,22:24) = S2(:,35:37);

if displays==1
    display('calibrated ENE changes');
    display(delta_eff(3,1:end));
end

output_covid(:,1) = output_covid(:,1) * 3600; % EUE $/KJ -> $/kWh
output_covid(:,2) = output_covid(:,2) / 3600; % EPE PJ/(t$)^0.3 / (billion cap)^0.7 -> PWh/(t$)^0.3 / (billion cap)^0.7
output_covid(:,12) = output_covid(:,12) / 3600; % energy PJ -> PWh
output_covid(:,21) = output_covid(:,21) / 3600; % cumulative green energy PJ -> PWh

end
