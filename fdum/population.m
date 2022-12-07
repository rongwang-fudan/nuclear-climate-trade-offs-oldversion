% Author: Rong Wang // contact rongwang@fudan.edu.cn //
% Date: 2021.2.10

function [L, Lcn]  =  population( LA )

global inputs realtime
%   realtime:  time
%   inputs 45x6: 1 energy PWh; 2 capital trill $; 3 GDP trill $; 4 population mill; 5 energy price ($/kWh); 6 omega

% Time horizon
T = size(realtime,1);

%Initial population (millions) 2015
% L0 = 7403;

%Population (millions)
L = zeros(T,1);

%Growth of population in 2015
x=realtime(41:45,1); y=inputs(41:45,4); [r,m,b] = regression(x',y');

%calibrated to meet the rate of growth in 2015 using data observed by IEA
Lg0 = log((inputs(45,4)+m)/inputs(45,4)) / log(LA/inputs(45,4));

for i=1:T % DICE-2013R
    if i<=size(inputs,1)
        L(i,1) = inputs(i,4) ;
    else
        L(i,1) = L(i-1,1) * (LA/L(i-1,1))^(Lg0 * realtime(i,2)) ;
    end
end

end
