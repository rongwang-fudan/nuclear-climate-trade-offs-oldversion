% Author: Rong Wang // contact rongwang@fudan.edu.cn //
% Date: 2022.12.6
tic
clear;

load('..\nuclear\nv.dat','-mat');
load('..\nuclear\scc.dat','-mat');

ncp=[0:0.1:1.5]; ncp(12:16)=[1.5 2 3 4 5]; ntw=ncp';
ntw2=[0.5 1 1.5 2 3 4]; ntw2=ntw2';
colorrb=[210 55 55
   255 88 137
   255 140 195
   255 188 225
   255 223 248
   234 215 255
   238 243 248
   212 212 255
   186 186 255
   156 156 255
   138 138 255
   83 83 255
   255 255 255]./255;
 colorrbline(:,1)=max(0,colorrb(:,1)-0.1);
 colorrbline(:,2)=max(0,colorrb(:,2)-0.1);
 colorrbline(:,3)=max(0,colorrb(:,3)-0.1);
 colorrbline(13,:)=[0 8 110]./255;
 
for exp=1:4
for i=13:-1:1
subplot(2,4,exp);
    plot(ntw(:,1),nv(:,(exp-1)*14+1+i),'-','LineWidth',1,'Color',colorrbline(14-i,1:3)); hold on;
subplot(2,4,exp+4);
    plot(ntw(2:15,1),scc(2:15,(exp-1)*14+1+i),'-','LineWidth',1,'Color',colorrbline(14-i,1:3)); hold on;

end
subplot(2,4,exp);
plot(ntw(:,1),nv(:,(exp-1)*14+1),'-','LineWidth',3,'Color',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','none','MarkerSize',6); hold on;
axis([0 5 -1500 1600]);

subplot(2,4,exp+4);
plot(ntw(2:15,1),scc(2:15,(exp-1)*14+1),'-','LineWidth',3,'Color',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','none','MarkerSize',6); hold on;
axis([0 4 -2000 1500]);
end