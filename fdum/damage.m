% Author: Rong Wang // contact rongwang@fudan.edu.cn //
% Date: 2021.6.18

function [ dcoef, xy_damage ] = damage( dpower, damagedata )

if damagedata==1
    xy = load('files\climate damage high.txt');
else
    % 126-160 2017-A Survey of Global Impacts of Climate Change: Replication, Survey Methods, and a Statistical Analysis
    xy = load('files\climate damage.txt');
end
% plot(xy(:,1),xy(:,2),'o');

fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0],...
               'Upper',[0.1],...
               'StartPoint',[0]);
if dpower==2
    ft =fittype('a*T1^2', 'independent', 'T1','options',fo);
elseif dpower==3
    ft =fittype('a*T1^3', 'independent', 'T1','options',fo);
end

[ps2,gof] = fit(xy(:,1),xy(:,2),ft);

n=size(xy,1);
xy_damage=zeros(59+n,3);
xy_damage(1:n,1)=xy(1:n,1);
xy_damage(1:n,3)=xy(1:n,2);

xy_damage((n+1):(n+59),1)=[0.2:0.1:6];
xy_damage(:,2)=ps2(xy_damage(:,1));

% plot(ps2, xy(:,1),xy(:,2));

dcoef = ps2(1);

display('Damage function: dpower for power coefficient on temperature, dcoef for damage as a percentage of GDP for 1 degree warming');
display([dcoef dpower]);

end

