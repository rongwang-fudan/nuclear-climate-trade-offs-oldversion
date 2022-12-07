% Author: Rong Wang
% Date: 2017.3.10

tic
clear;

ene = load('files\IEA_TPES.txt'); % 161x44, Mtoe, from 1971 to 2014
ene = ene.*(66464/1587.5); % Mtoe to PJ
pop = load('files\IEA_POP.txt'); % million
gdp = load('files\IEA_GDPppp.txt'); % billion 2005 $, DICE-R IEA_GDPexc.txt IEA_GDPppp
capital = load('files\capital_stock_1970_2014.txt'); % 46 *5 china, usa, eu, japan, germany; 1970 - 2015
us_eneshare = load('files\us_energyshare2015.txt'); % 46x1 for 1970 to 2015
sweden_pe = load('files\sweden_pe.txt');
newdata=[2008.113	2017.424	32714.676	33644.921	7324.4	7354.3	7413
1369.435	1376.049	5390.257	5762.185	2970.6	3005.9	3053
319.449	321.774	14701.679	15083.356	2296.5	2275.9	2272.7
192.434	193.028	7664.69	7780.87	1605	1626.7	1642
126.795	126.573	4958.05	5018.51	452.3	445.8	445.3
1295.3	1311.05	1559.396	1677.339	663.6	685.1	723.9
]; % pop14 pop15; dgp14 dgp 15; ene14 ene15  -- globe china us eu japan india
cninfo = load('files\country_info.txt');
input=zeros(45,31);
pesce = 1; % 1 for us price; 2 for sweden price
if pesce==1
    for j=1:44
        me = us_eneshare(j+1) / 100 * gdp(4,j) / ene(4,j); % energy expenditures as share of GDP in the US from 1971 to 2014
        input(j,1) = me*3.6; % $/kWh
    end
    input(45,1) = us_eneshare(46) / 100 * gdp(4,44) * newdata(3,4) / newdata(3,3) / (ene(4,44) * newdata(3,6) / newdata(3,5) )*3.6; % $/kWh
else
    input(1:45,1) = sweden_pe(1:45,1);
end
for j=1:44
    e0=0; y0=0; p0=0;
    me=input(j,1)/3.6; % $/kWh to $/MJ
    for i=1:162
        if i<=161
            if ene(i,j)<0 || pop(i,j)<0 || gdp(i,j)<0
                continue;
            end
            if cninfo(i)<10
                e0=e0+ene(i,j); y0=y0+gdp(i,j); p0=p0+pop(i,j);
            end
            if i==120
                cn=1; % china
            elseif i==4
                cn=2; % usa
            elseif i==37
                cn=3; % eu
            elseif i==8
                cn=4; % japan
            elseif i==103
                cn=5; % India
            else
                cn=0;
            end
            if cn~=0
                input(j,cn+2) = ene(i,j)/pop(i,j)/3.6;
                input(j,cn+8) = gdp(i,j)/pop(i,j);
                input(j,cn+14) = ene(i,j)/gdp(i,j)*me;
                input(j,cn+26) = pop(i,j);
            end
        else
            input(j,2) = e0/p0/3.6;
            input(j,8) = y0/p0;
            input(j,14) = e0*me/y0;
            input(j,26) = p0;
        end
    end
end
input(45,26) = input(44,26)*2-input(43,26);

for i=1:6
    input(45,i+1) = input(44,i+1) * newdata(i,6) / newdata(i,5) / (newdata(i,2) / newdata(i,1));
    input(45,i+7) = input(44,i+7) * newdata(i,4) / newdata(i,3) / (newdata(i,2) / newdata(i,1));
    input(45,i+13) = input(45,i+1) * input(45,1) / input(45,i+7);
end

for i=1:44
    input(i,20) = (capital(i,1)+capital(i,2)+capital(i,3)+capital(i,4))/(pop(120,i)+pop(4,i)+pop(37,i)+pop(8,i));
    input(i,21) = capital(i,1)/pop(120,i);
    input(i,22) = capital(i,2)/pop(4,i);
    input(i,23) = capital(i,3)/pop(37,i);
    input(i,24) = capital(i,4)/pop(8,i);
end
popo1=pop(120,44)*newdata(2,2)/newdata(2,1); input(45,21) = capital(45,1)/popo1;
popo2=pop(4,44)*newdata(3,2)/newdata(3,1); input(45,22) = capital(45,2)/popo2;
popo3=pop(37,44)*newdata(4,2)/newdata(4,1); input(45,23) = capital(45,3)/popo3;
popo4=pop(8,44)*newdata(5,2)/newdata(5,1); input(45,24) = capital(45,4)/popo4;
input(45,20) = (capital(45,1)+capital(45,2)+capital(45,3)+capital(45,4))/(popo1+popo2+popo3+popo4);

% china us eu japan india
input(45,27) = popo1;
input(45,28) = popo2;
input(45,29) = popo3;
input(45,30) = popo4;
input(45,31) = input(44,31)*2-input(43,31);

n=size(input,2);
for i=1:n
    input(:,i)=smooth(input(:,i),'sgolay',2);
end


lra=zeros(26*5,6);
for id=1:5
    
for i=1:26
    if i<=26
        x=[(1970+i):(1989+i)]; y=input(i:(i+19),1); [sR,lr_pe0,bb0] = regression(x,log(y'));
        y=input(i:(i+19),1+id); [sR,lr_e0,bb0] = regression(x,log(y'));
        y=input(i:(i+19),7+id); [sR,lr_y0,bb0] = regression(x,log(y'));
        y=input(i:(i+19),13+id); avef0=mean(y,1); startf0=input(i,13+id); startpe0=input(i,1);
        y=input(i:(i+19),19+id); [sR,lr_k0,bb0] = regression(x,log(y'));
    else
        x=[(1970+i):2015]; y=input(i:45,1); [sR,lr_pe0,bb0] = regression(x,log(y'));
        y=input(i:45,1+id); [sR,lr_e0,bb0] = regression(x,log(y'));
        y=input(i:45,7+id); [sR,lr_y0,bb0] = regression(x,log(y'));
        y=input(i:45,13+id); avef0=mean(y,1); startf0=input(i,13+id); startpe0=input(i,1);
        y=input(i:45,19+id); [sR,lr_k0,bb0] = regression(x,log(y'));
    end
    lr_se0 = lr_pe0 + lr_e0 - lr_y0;
    
    lra(i+(id-1)*26,1) = lr_se0;
    lra(i+(id-1)*26,2) = lr_e0;
    lra(i+(id-1)*26,3) = lr_k0;
    lra(i+(id-1)*26,4) = lr_y0;
    lra(i+(id-1)*26,5) = avef0;
    lra(i+(id-1)*26,6) = startf0;

end

end

save('files\EconomicHistoricalData_1971_2015.dat','input');
% 45x31; 1 energy price ($/kWh)
% 2-7 energy per capita (MWh/cap) for globe china us eu japan india
% 8-13 gdp per capita (MWh/cap) for globe china us eu japan india
% 14-19 omega for globe china us eu japan india
% 20-25 capital stock per capita (k$/cap) for globe china us eu japan india
% 26-31 population (million) for globe china us eu japan india



