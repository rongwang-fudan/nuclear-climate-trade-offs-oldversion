% Author: Rong Wang // contact rongwang@fudan.edu.cn //
% Date: 2021.6.18
% Calibtration of calrsav by capital, energy and ouput
%    savings rate 0.25 for default

function [ oilprice1971_2019, oilprice2020_normalize, energyprice2020, lostfunction_SE, scatterplots, co2emi_2019_2021] = energyprice( covidyear )

global inputs
%   inputs 45x6: 1 energy PWh; 2 capital trill $; 3 GDP trill $; 4 population mill; 5 energy price ($/kWh); 6 omega

inp = load('files\energyprices.txt');
% monthly energy prices from different sources
% 1	Coal, Australian thermal coal - Monthly Price		
% 2	Coal, South African export price - Monthly Price		
% 3	Coal, Colombia - Monthly Price		
% 4	Crude Oil (petroleum) - Monthly Price		
% 5	Crude Oil (petroleum); Dated Brent - Monthly Price		
% 6	Crude Oil (petroleum); Dubai Fateh - Monthly Price		
% 7	Crude Oil (petroleum); West Texas Intermediate - Monthly Price		
% 8	Diesel - Monthly Price		
% 9	Gasoline - Monthly Price		
% 10	Heating Oil - Monthly Price		
% 11	Indonesian Liquified Natural Gas - Monthly Price		
% 12	Jet Fuel - Monthly Price		
% 13	Natural Gas - Monthly Price		
% 14	Propane - Monthly Price		
% 15	RBOB Gasoline - Monthly Price		
% 16	Russian Natural Gas - Monthly Price
oilprice2020 = zeros(29,9);
es = zeros(49,18);
ids = [4 5 6 7 8 9 10 12 15];

for k=1:9
    for i=1:size(inp,1)
        if inp(i,1)>=(covidyear-1) && inp(i,4)==ids(k)
            oilprice2020(inp(i,2)+(inp(i,1)-2019)*12,k)=inp(i,3);
        end
        if inp(i,4)==ids(k) && inp(i,1)<=2019
            es(inp(i,1)-1970,k)=es(inp(i,1)-1970,k)+inp(i,3);
            es(inp(i,1)-1970,k+9)=es(inp(i,1)-1970,k+9)+1;
        end
    end
end

oilprice2020_normalize = zeros(29,9);
for i=1:9
    ref=mean(oilprice2020(1:12,i),1);
    oilprice2020_normalize(:,i)=oilprice2020(:,i)./ref;
end

oilprice1971_2019 = zeros(49,9);
for i=1:49
    for j=1:9
        if es(i,j+9)>0
            oilprice1971_2019(i,j)=es(i,j)/es(i,j+9);
        else
            oilprice1971_2019(i,j)=-999;
        end
    end
end

lostfunction_SE = zeros(36,10);
scatterplots = zeros(198,37);
for i=1:36
    e0=inputs(:,5)-(i-1)*0.001; % energy price minus the fixed cost
    n=0;
    n1=1;
    for k=1:9
        idx=find(oilprice1971_2019(:,k)>0);
        % market oil price
        p1=oilprice1971_2019(idx(1):45,k);
        p2=p1./mean(p1,1); % the average normalized to 1
        % energy price
        e2=e0(idx(1):45)./mean(e0(idx(1):45),1); % the average normalized to 1
        % calculate the mean squared deviation between the two prices
        lostfunction_SE(i,k)=sum((e2-p2).*(e2-p2),1)./(46-idx(1));
        lostfunction_SE(i,10)=lostfunction_SE(i,10)+sum((e2-p2).*(e2-p2),1);
        n=n+46-idx(1);
        % list the original data out for a scatter plot in Excel
        n2=n1+45-idx(1);
        scatterplots(n1:n2,1)=p2;
        scatterplots(n1:n2,i+1)=e2;
        n1=n2+1;
    end
    lostfunction_SE(i,10)=lostfunction_SE(i,10)/n;
end

energyprice2020=zeros(29,9);
for k=1:9
    [B,IX] = sort(lostfunction_SE(:,k),1);
    fixcost = (IX(1)-1)*0.001;
    for i=1:29
        energyprice2020(i,k) = (inputs(45,5)-fixcost)/oilprice1971_2019(45,k)*oilprice2020(i,k) + fixcost;
    end
end

co2emi = load('files\co2emi_2019_2021.txt');
% 1 country; 2 sector; 3 year; 4 month; 5 emission MtCO2/day
% sector: 1	Power, 2 Industry, 3 Ground Transport, 4 Residential, 5 Dom Aviation, 6 Inter Aviation, 7 Inter Shipping, 8 Total
% country: 1 Brazil, 2 China,3 EU27 & UK, 4 France,5,Germany, 6 India, 7 Italy, 8 Japan, 9 ROW, 10 Russia, 11 Spain, 12 UK, 13 US, 14 WORLD
co2emi_2019_2021 = zeros(29,8);
for i=1:size(co2emi,1)
    if co2emi(i,3)>=2019
        mm=floor(co2emi(i,4)+(co2emi(i,3)-2019)*12);
        if co2emi(i,1)==14
            co2emi_2019_2021(mm,co2emi(i,2))=co2emi_2019_2021(mm,co2emi(i,2))+co2emi(i,5)*12/1e3; % GtCO2/year
        end
    end
end

end