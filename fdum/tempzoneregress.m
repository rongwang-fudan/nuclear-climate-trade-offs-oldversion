% Author: Rong Wang // contact rongwang@fudan.edu.cn //
% Date: 2021.2.10

function regress_tempzone  =  tempzoneregress( regionids , plotreg )

% regression of regional temperature with global temperature
regress_tempzone=zeros(regionids,4); %1  N Ame; 2 S Ame; 3 W Eur; 4 N Afr & Middle East; 5 S Afr; 6 E Eur; 7 E Asia; 8 S & SE Asia; 9 Oceania; 10 Polar; 11 Ocean
tempinput = load('files\temperature.txt'); % 764x12 (1 globe; 2-11 region; 12 for ocean)
tempinputyy=zeros(63,12);
for yy=1:63
    tempinputyy(yy,1:12)=mean(tempinput((yy*12-11):(yy*12),1:12),1);
end
for regionid=1:regionids
    [b2,bint2,r2,rint2,stats2]=regress(tempinputyy(:,regionid+1),[ones(63,1) tempinputyy(:,1)]);
    regress_tempzone(regionid,1)=b2(2,1); % slope
    regress_tempzone(regionid,2)=(bint2(2,2)-bint2(2,1))/1.96/2; % std of slope
    regress_tempzone(regionid,3)=mean(tempinputyy(:,regionid+1),1); % ave of y
    regress_tempzone(regionid,4)=mean(tempinputyy(:,1),1); % ave of x
    
    reg_statis(regionid,1)=b2(2,1); % slope
    reg_statis(regionid,2)=regress_tempzone(regionid,3)-b2(2,1)*regress_tempzone(regionid,4);
    reg_statis(regionid,3)=stats2(1); % R2
    reg_statis(regionid,4)=stats2(3); % P

end

if plotreg==1
for regp=2:12
    subplot(3,4,regp-1);
    y=tempinputyy(:,regp);
    x=tempinputyy(:,1);
    a=polyfit(x,y,1); 
    x2=(min(tempinputyy(:,1)):0.1:max(tempinputyy(:,1)));
    y2=polyval(a,x2);
    
    [p,S]=polyfit(x,y,1);
    alpha = 0.05; % Significance level
   [yfit,delta] = polyconf(p,x2,S,'alpha',alpha);
   y3=yfit-delta;
   y4=yfit+delta;
%     plot(x2,yfit-delta,'LineStyle','-','LineWidth',1,'Color','r'); hold on; 
%     plot(x2,yfit+delta,'LineStyle','-','LineWidth',1,'Color','r'); hold on; 
     pic01=fill([x2,fliplr(x2)],[y3,fliplr(y4)],[0.8 0.8 0.8]); hold on;
     set(pic01,'edgealpha', 0, 'facealpha', 0.4);
    plot(tempinputyy(:,1),tempinputyy(:,regp),'.','Color',[0 0 0],'MarkerSize',8);hold on;
    plot(x2,y2,'LineStyle','-','LineWidth',2,'Color','b'); hold on; 

end

end

end
