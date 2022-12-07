% Author: Rong Wang // contact rongwang@fudan.edu.cn //
% Date: 2021.6.18

function [iec_cn] = Calibration_IEC( cndata, plotiec )

% Economic variables cndata(KEYL)
% plotiec=1 make a plot
% Economic variables econcn=zeros(49,23,i2+1)
% 1 EUE; 2 EPE; 3 ENE; 4 backstop price $/tCO2
% 5 abatement cost as a percentage of GDP; 6 climate damage as a percentage of GDP; 7 net output
% 8 fraction of labor allocation to energy; 9 fraction of investment allocation to energy
% 10 total capital (t$); 11 energy capital (t$); 12 energy (PJ); 13 output (t$); 14 energy price $/kWh; 15 omega
% 16 fraction of energy investment allocated to zero-carbon energy
% 17 energy capital carbon-emission-free (t$);
% 18 fraction to abate CO2 emission; 19 carbon price $/tCO2; 20 CO2 emissions Gt CO2; 21 green energy PJ; 22 invest change, 23 labor change

global alpha elas inputs cou_iform
%   alpha:  elasticity of output to capital
%   elas:   elasticity of substitution between energy and non-energy in output

cn_num=size(cndata,1);
% piecewise=1; % 1 use piecewise regression; 2 use simple regression

% 1. Eastern and Southern Africa; 2. Northern Africa; 3. Western and Central Africa; 4. East Asia; 5. South and South-east Asia; 6. Western and Central Asia;
% 7. Europe; 8. Caribbean; 9. Central America; 10. North America; 11. Oceania; 12. South America
zonenum=8;
z8data=zeros(zonenum,49,4); % 1 Globe; 2 Africa, 3 East Asia; 4 South/Southeast Asia; 5 Middle East; 6 Europe; 7 North America; 8 Oceania
z8id=[2 2 2 3 4 5 6 0 0 7 8 0];
for i=1:49
    for k=1:4
        z8data(1,i,k)=cndata(1,(k*49-48+i)); % Globe
        for cn=2:cn_num
            z8=z8id(1,cou_iform(cndata(cn,1),3));
            if z8>=2 && z8<=zonenum
                z8data(z8,i,k)=z8data(z8,i,k)+cndata(cn,(k*49-48+i));
            end
        end
    end
end

xy_iec_zone=zeros(30,5,zonenum);
iec_statis=zeros(zonenum,12);
slopezone=zeros(zonenum,28);
for cn=1:zonenum
lra=zeros(30,7);
for i=1:30
    x=[(1970+i):(1989+i)];
    x15=[(1970+i):min(1989+i,2015)];
    i15=min((i+19),45);
    y=inputs(i:i15,5); [sR,lr_pe0,bb0] = regression(x15,log(y')); % energy price
    y=z8data(cn,i:(i+19),2)./z8data(cn,i:(i+19),4); [sR,lr_e0,bb0] = regression(x,log(y)); % rate of e=E/L
    y=z8data(cn,i:(i+19),3)./z8data(cn,i:(i+19),4); [sR,lr_y0,bb0] = regression(x,log(y)); % rate of y=Y/L
    y=z8data(cn,i:(i+19),1)./z8data(cn,i:(i+19),4); [sR,lr_k0,bb0] = regression(x,log(y)); % rate of k=K/L    
    y=zeros(i15-i+1,1);
    for j=i:i15
        y(j-i+1,1)=inputs(j,6) * (z8data(cn,j,2)/z8data(cn,j,3)) / (inputs(j,1)*3600/inputs(j,3)); % omega    
    end
    avef0=mean(y,1); startf0=y(1);
    lr_se0 = lr_pe0 + lr_e0 - lr_y0;
    lr_B = lr_e0 - alpha*lr_k0;
    lr_A = lr_y0 - alpha*lr_k0;
    lr_taue0 = elas/(elas-1)*lr_se0 - lr_e0 + lr_y0;
    lr_we0 = lr_B - lr_se0;
    lr_b0 = (lr_A - avef0*(lr_we0 + lr_taue0))/(1-avef0);
    lra(i,1) = startf0;
    lra(i,2) = lr_taue0;
    lra(i,3) = lr_we0;
    lra(i,4) = max(0.0001,lr_b0);
    lra(i,5) = 0; % damage of climate change
    lra(i,6) = lra(i,1);
    lra(i,7) = 1-lra(i,1);
end
xy_iec_zone(:,1,cn)=lra(:,2); % eue rate
xy_iec_zone(:,2,cn)=lra(:,3); % epe rate
xy_iec_zone(:,3,cn)=lra(:,4); % ene rate
xy_iec_zone(:,4,cn)=lra(:,6); % omega
xy_iec_zone(:,5,cn)=lra(:,7); % 1 - omega

% EUE
x=log10(lra(:,6)); y=lra(:,2);
% if piecewise==1
%     slopezone(cn,1) = (max(x,[],1)+min(x,[],1))/2; % divide the data into two parts
% else
    slopezone(cn,1) = 0;
% end
idx1=find(x<slopezone(cn,1));
[b2,bint2,r2,rint2,stats2]=regress(y(idx1,1),[ones(size(idx1,1),1) x(idx1,1)]);
slopezone(cn,2) = b2(2);
slopezone(cn,3) = abs((bint2(2,2)-bint2(2,1))/1.96/2/b2(2));
slopezone(cn,16) = mean(x(idx1,1),1);
slopezone(cn,17) = mean(y(idx1,1),1);
iec_statis(cn,1)=stats2(1); % R2
iec_statis(cn,2)=stats2(3); % P
if slopezone(cn,1)<-0.001
    idx1=find(x>=slopezone(cn,1));
    [b2,bint2,r2,rint2,stats2]=regress(y(idx1,1),[ones(size(idx1,1),1) x(idx1,1)]);
    slopezone(cn,4) = b2(2);
    slopezone(cn,5) = abs((bint2(2,2)-bint2(2,1))/1.96/2/b2(2));    
    slopezone(cn,18) = mean(x(idx1,1),1);
    slopezone(cn,19) = mean(y(idx1,1),1);
    iec_statis(cn,3)=stats2(1); % R2
    iec_statis(cn,4)=stats2(3); % P
end

% EPE
y=lra(:,3);
% if piecewise==1
%     slopezone(cn,6) = (max(x,[],1)+min(x,[],1))/2; % divide the data into two parts
% else
    slopezone(cn,6) = 0;
% end
idx1=find(x<slopezone(cn,6));
[b2,bint2,r2,rint2,stats2]=regress(y(idx1,1),[ones(size(idx1,1),1) x(idx1,1)]);
slopezone(cn,7) = b2(2);
slopezone(cn,8) = abs((bint2(2,2)-bint2(2,1))/1.96/2/b2(2));
slopezone(cn,20) = mean(x(idx1,1),1);
slopezone(cn,21) = mean(y(idx1,1),1);
iec_statis(cn,5)=stats2(1); % R2
iec_statis(cn,6)=stats2(3); % P
if slopezone(cn,6)<-0.001
    idx1=find(x>=slopezone(cn,6));
    [b2,bint2,r2,rint2,stats2]=regress(y(idx1,1),[ones(size(idx1,1),1) x(idx1,1)]);
    slopezone(cn,9) = b2(2);
    slopezone(cn,10) = abs((bint2(2,2)-bint2(2,1))/1.96/2/b2(2));    
    slopezone(cn,22) = mean(x(idx1,1),1);
    slopezone(cn,23) = mean(y(idx1,1),1);
    iec_statis(cn,7)=stats2(1); % R2
    iec_statis(cn,8)=stats2(3); % P
end

% ENE - piecewise
idxene=find(lra(:,4)>=0.0015); % disturbed economy
x=log10(lra(idxene,7));
y=lra(idxene,4);
yrange=max(y,[],1)-min(y,[],1);
if yrange>0.01
    slopezone(cn,11) = (max(x,[],1)+min(x,[],1))/2; % divide the data into two parts
    slopezone(cn,11) = min(slopezone(cn,11), prctile(x,[max(0,100-floor(18/size(idxene,1)*100))])); % ensure 60% of 30 data points
else
    slopezone(cn,11) = 0;
end
idx1=find(x<slopezone(cn,11));
if size(idx1)<3
    slopezone(cn,11) = 0;
end
idx1=find(x<slopezone(cn,11));
[b2,bint2,r2,rint2,stats2]=regress(y(idx1,1),[ones(size(idx1,1),1) x(idx1,1)]);
slopezone(cn,12) = max(0.1,b2(2));
slopezone(cn,13) = abs((bint2(2,2)-bint2(2,1))/1.96/2/b2(2));
slopezone(cn,24) = mean(x(idx1,1),1);
slopezone(cn,25) = mean(y(idx1,1),1);
iec_statis(cn,9) = stats2(1);  % R2
iec_statis(cn,10) = stats2(3); % P
if slopezone(cn,11)<-0.001
    idx1=find(x>=slopezone(cn,11));
    [b2,bint2,r2,rint2,stats2]=regress(y(idx1,1),[ones(size(idx1,1),1) x(idx1,1)]);
    slopezone(cn,14) = max(0.1,b2(2));
    slopezone(cn,15) = abs((bint2(2,2)-bint2(2,1))/1.96/2/b2(2));    
    slopezone(cn,26) = mean(x(idx1,1),1);
    slopezone(cn,27) = mean(y(idx1,1),1);
    iec_statis(cn,11)=stats2(1); % R2
    iec_statis(cn,12)=stats2(3); % P
end
% Use a simple regression if the piecewise regression fails
if slopezone(cn,11)<-0.001 && iec_statis(cn,12)>0.05
    slopezone(cn,11) = 0;
    idx1=find(x<slopezone(cn,11));
    [b2,bint2,r2,rint2,stats2]=regress(y(idx1,1),[ones(size(idx1,1),1) x(idx1,1)]);
    slopezone(cn,12) = max(0.1,b2(2));
    slopezone(cn,13) = abs((bint2(2,2)-bint2(2,1))/1.96/2/b2(2));
    slopezone(cn,14:15) = 0;
    slopezone(cn,24) = mean(x(idx1,1),1);
    slopezone(cn,25) = mean(y(idx1,1),1);
    iec_statis(cn,9)=stats2(1);  % R2
    iec_statis(cn,10)=stats2(3); % P
end

% Determine the latest omega
slopezone(cn,28)=lra(end,6);

end

if plotiec==1
    % EUE rates as a function of log10(Omega)
    for cn=1:zonenum
        subplot(6,4,cn);
        x0=log10(xy_iec_zone(:,4,cn)); 
        idx0=find(x0<slopezone(cn,1));
        x=log10(xy_iec_zone(idx0,4,cn)); y=xy_iec_zone(idx0,1,cn);
        a=polyfit(x,y,1); 
        x2=[-1.7:0.1:-0.8]; y2=polyval(a,x2);
        i2=size(x2,2); y3=zeros(2,i2);
        for i=1:i2
            y3(1,i) = mean(y,1) + (y2(1,i)-mean(y,1)) * (1-slopezone(cn,3)*1.96);
            y3(2,i) = mean(y,1) + (y2(1,i)-mean(y,1)) * (1+slopezone(cn,3)*1.96);
        end
        plot(x,y,'o','MarkerEdgeColor',[0 0.6 0.6],'MarkerFaceColor','none','MarkerSize',3); hold on;   
        plot(x2,y2,'LineStyle','-','LineWidth',2,'Color',[0 0.6 0.6]); hold on;        
        plot(x2,y3(1,:),'LineStyle',':','LineWidth',1,'Color',[0 0.6 0.6]); hold on;        
        plot(x2,y3(2,:),'LineStyle',':','LineWidth',1,'Color',[0 0.6 0.6]); hold on;
        %
        if slopezone(cn,1)<-0.001
        idx0=find(x0>=slopezone(cn,1));
        x=log10(xy_iec_zone(idx0,4,cn)); y=xy_iec_zone(idx0,1,cn);
        a=polyfit(x,y,1); 
        x2=[-1.7:0.1:-0.8]; y2=polyval(a,x2);
        i2=size(x2,2); y3=zeros(2,i2);
        for i=1:i2
            y3(1,i) = mean(y,1) + (y2(1,i)-mean(y,1)) * (1-slopezone(cn,5)*1.96);
            y3(2,i) = mean(y,1) + (y2(1,i)-mean(y,1)) * (1+slopezone(cn,5)*1.96);
        end
        plot(x,y,'o','MarkerEdgeColor',[1 0 0],'MarkerFaceColor','none','MarkerSize',3); hold on;        
        plot(x2,y2,'LineStyle','-','LineWidth',2,'Color',[0.8 0 0]); hold on;        
        plot(x2,y3(1,:),'LineStyle',':','LineWidth',1,'Color',[0.8 0 0]); hold on;        
        plot(x2,y3(2,:),'LineStyle',':','LineWidth',1,'Color',[0.8 0 0]); hold on;
        end
        axis([-1.7 -0.8 -0.08 0.08]);
    end
    
    % EPE rates as a function of log10(Omega)
    for cn=1:zonenum
        subplot(6,4,cn+8);
        x0=log10(xy_iec_zone(:,4,cn)); 
        idx0=find(x0<slopezone(cn,6));
        x=log10(xy_iec_zone(idx0,4,cn)); y=xy_iec_zone(idx0,2,cn);
        a=polyfit(x,y,1); 
        x2=[-1.7:0.1:-0.8]; y2=polyval(a,x2);
        i2=size(x2,2); y3=zeros(2,i2);
        for i=1:i2
            y3(1,i) = mean(y,1) + (y2(1,i)-mean(y,1)) * (1-slopezone(cn,8)*1.96);
            y3(2,i) = mean(y,1) + (y2(1,i)-mean(y,1)) * (1+slopezone(cn,8)*1.96);
        end
        plot(x,y,'o','MarkerEdgeColor',[0.6 0.2 1],'MarkerFaceColor','none','MarkerSize',3); hold on;   
        plot(x2,y2,'LineStyle','-','LineWidth',2,'Color',[0.6 0.2 1]); hold on;        
        plot(x2,y3(1,:),'LineStyle',':','LineWidth',1,'Color',[0.6 0.2 1]); hold on;        
        plot(x2,y3(2,:),'LineStyle',':','LineWidth',1,'Color',[0.6 0.2 1]); hold on;
        %
        if slopezone(cn,6)<-0.001
        idx0=find(x0>=slopezone(cn,6));
        x=log10(xy_iec_zone(idx0,4,cn)); y=xy_iec_zone(idx0,2,cn);
        a=polyfit(x,y,1); 
        x2=[-1.7:0.1:-0.8]; y2=polyval(a,x2);
        i2=size(x2,2); y3=zeros(2,i2);
        for i=1:i2
            y3(1,i) = mean(y,1) + (y2(1,i)-mean(y,1)) * (1-slopezone(cn,10)*1.96);
            y3(2,i) = mean(y,1) + (y2(1,i)-mean(y,1)) * (1+slopezone(cn,10)*1.96);
        end
        plot(x,y,'o','MarkerEdgeColor',[1 0 0],'MarkerFaceColor','none','MarkerSize',3); hold on;        
        plot(x2,y2,'LineStyle','-','LineWidth',2,'Color',[0.8 0 0]); hold on;        
        plot(x2,y3(1,:),'LineStyle','--','LineWidth',1,'Color',[1 0.7 0.7]); hold on;        
        plot(x2,y3(2,:),'LineStyle','--','LineWidth',1,'Color',[1 0.7 0.7]); hold on;
        end
        axis([-1.7 -0.8 -0.08 0.08]);
    end
    
    % ENE rates as a function of log10(1-Omega)
    for cn=1:zonenum
        subplot(6,4,cn+16);
        idxene=find(xy_iec_zone(:,3,cn)>=0.0015); % destroyed economy
        x0=log10(xy_iec_zone(idxene,5,cn)); y0=xy_iec_zone(idxene,3,cn);
        idx0=find(x0<slopezone(cn,11));
        x=x0(idx0,1); y=y0(idx0,1);
        a=polyfit(x,y,1); 
        x2=[-0.08:0.01:0]; y2=polyval(a,x2);
        i2=size(x2,2); y3=zeros(2,i2);
        for i=1:i2
            y3(1,i) = mean(y,1) + (y2(1,i)-mean(y,1)) * (1-slopezone(cn,13)*1.96);
            y3(2,i) = mean(y,1) + (y2(1,i)-mean(y,1)) * (1+slopezone(cn,13)*1.96);
        end
        plot(x,y,'o','MarkerEdgeColor',[0 0 1],'MarkerFaceColor','none','MarkerSize',3); hold on;
        if slopezone(cn,12)>0.1 || cn==7
            plot(x2,y2,'LineStyle','-','LineWidth',2,'Color',[0 0 1]); hold on;        
            plot(x2,y3(1,:),'LineStyle',':','LineWidth',1,'Color',[0.64 0.72 1]); hold on;        
            plot(x2,y3(2,:),'LineStyle',':','LineWidth',1,'Color',[0.64 0.72 1]); hold on;  
        end
        %
        if slopezone(cn,11)<-0.001
        idx0=find(x0>=slopezone(cn,11));
        x=x0(idx0,1); y=y0(idx0,1);
        a=polyfit(x,y,1); 
        x2=[-0.08:0.01:0]; y2=polyval(a,x2);
        i2=size(x2,2); y3=zeros(2,i2);
        for i=1:i2
            y3(1,i) = mean(y,1) + (y2(1,i)-mean(y,1)) * (1-slopezone(cn,15)*1.96);
            y3(2,i) = mean(y,1) + (y2(1,i)-mean(y,1)) * (1+slopezone(cn,15)*1.96);
        end
        plot(x,y,'o','MarkerEdgeColor',[0.8 0 0],'MarkerFaceColor','none','MarkerSize',3); hold on;     
        plot(x2,y2,'LineStyle','-','LineWidth',2,'Color',[0.8 0 0]); hold on;        
        plot(x2,y3(1,:),'LineStyle',':','LineWidth',1,'Color',[1 0.7 0.7]); hold on;        
        plot(x2,y3(2,:),'LineStyle',':','LineWidth',1,'Color',[1 0.7 0.7]); hold on;
        end
        % outliers
        idxene2=find(xy_iec_zone(:,3,cn)<0.0015); % disturbed economy
        plot(log10(xy_iec_zone(idxene2,5,cn)),xy_iec_zone(idxene2,3,cn),'o','MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor','none','MarkerSize',3); hold on;  
        axis([-0.066 -0.01 -0.002 0.052]);
    end
end

% Smooth the curve
for cn=1:zonenum
    if slopezone(cn,11)<-0.001
        slopezone(cn,11)=(slopezone(cn,25)-slopezone(cn,27)+slopezone(cn,26)*slopezone(cn,14)-slopezone(cn,24)*slopezone(cn,12))/(slopezone(cn,14)-slopezone(cn,12));
    else
        slopezone(cn,11)=-0.03;
        slopezone(cn,14)=max(0.1,slopezone(cn,12));
        slopezone(cn,15)=slopezone(cn,13);
        slopezone(cn,12)=0.1;
    end
end

% Data by country
iec_cn=zeros(cn_num,29);
for cn=1:cn_num
    if cn==1
        iec_cn(cn,1:28)=slopezone(1,1:28); % Globe
        continue;
    end
    z8=z8id(1,cou_iform(cndata(cn,1),3));
    if cndata(cn,1)==99 || cndata(cn,1)==105
        z8=7;
    end
    if cndata(cn,1)==129 || z8==0
        z8=3; % for mexico
    end
    if cou_iform(cndata(cn,1),3)==12
        z8=4; % for SA
    end
    iec_cn(cn,1:28)=slopezone(z8,1:28);
    
    % define the emerging countries
    if iec_cn(cn,14)>0.5 && cou_iform(cndata(cn,1),6)>14
        iec_cn(cn,29)=1;
    end
end

end


