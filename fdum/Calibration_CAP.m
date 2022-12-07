% Author: Rong Wang // contact rongwang@fudan.edu.cn //
% Date: 2021.6.18
% Finding the optimal calrsav to achieve the best agreement with observations
% Calibtration of calrsav 1 EUE; 2 ENE; 3-5 saving rate 1971-2003; 2004-2008; 2009-2015
%    savings rate 0.25 for default

function [ calrsav, output_cap ] = Calibration_CAP( L, iec, LR, switcher, displays, sensitivity )

global S0

calrsav0 = [0.9500    1.1000    0.2000    0.3500    0.3100    0.3100];
[S2, LFmin] = capital( L, iec, LR, switcher, calrsav0 );
calrsav2 = calrsav0;

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

for z=1:5
    
calrsavnew=calrsav0;
calrsavnew(z)=calrsavnew(z)*1.1;
[ss, LF] = capital( L, iec, LR, switcher, calrsavnew );
if LF(1)<=LFmin(1) && LF(2)<=LFmin(2) && LF(3)<=LFmin(3) && LF(4)<=LFmin(4)
    while LF(1)<=LFmin(1) && LF(2)<=LFmin(2) && LF(3)<=LFmin(3) && LF(4)<=LFmin(4) && (calrsav0(z)/calrsav2(z))<2
        if displays==1
            display('z+1');
        end
        S2=ss;
        calrsav0=calrsavnew;
        LFmin=LF;
        % try to use new parameters
        calrsavnew(z)=calrsavnew(z)*1.1;
        [ss, LF] = capital( L, iec, LR, switcher, calrsavnew );
    end
else
    calrsavnew=calrsav0;
    calrsavnew(z)=calrsavnew(z)*0.9;
    [ss, LF] = capital( L, iec, LR, switcher, calrsavnew );
    while LF(1)<=LFmin(1) && LF(2)<=LFmin(2) && LF(3)<=LFmin(3) && LF(4)<=LFmin(4) && (calrsav0(z)/calrsav2(z))>0.5
        if displays==1
            display('z-1');
        end
        S2=ss;
        calrsav0=calrsavnew;
        LFmin=LF;
        % try to use new parameters
        calrsavnew(z)=calrsavnew(z)*0.9;
        [ss, LF] = capital( L, iec, LR, switcher, calrsavnew );
    end
end

end

calrsav = calrsav0;
S0(1:45,1:23) = S2(1:45,1:23);

output_cap = S2;

if displays==1
    display('calibrated saving rate - calrsav: 1 for 1971-2003; 2009-2015; 2 for 2004-2008; 3 for 2009-2015');
    display(calrsav(3:5));
end

output_cap(:,1)  = output_cap(:,1) * 3600; % EUE $/KJ -> $/kWh
output_cap(:,24) = output_cap(:,24) * 3600; % EUE $/KJ -> $/kWh
output_cap(:,2)  = output_cap(:,2) / 3600; % EPE PJ/(t$)^0.3 / (billion cap)^0.7 -> PWh/(t$)^0.3 / (billion cap)^0.7
output_cap(:,25) = output_cap(:,25) / 3600; % EPE PJ/(t$)^0.3 / (billion cap)^0.7 -> PWh/(t$)^0.3 / (billion cap)^0.7
output_cap(:,12) = output_cap(:,12) / 3600; % energy PJ -> PWh
output_cap(:,21) = output_cap(:,21) / 3600; % cumulative green energy PJ -> PWh
output_cap(:,28) = output_cap(:,28) / 3600; % energy PJ -> PWh

end
