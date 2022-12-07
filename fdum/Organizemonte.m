% Author: Rong Wang // contact rongwang@fudan.edu.cn //
% Date: 2022.12.6
tic
clear;

load('files\cndata134.dat','-mat'); cndata=cndata2; clear cndata2;
load('..\nuclear\output_pop.dat','-mat'); % output_pop=zeros(400,cn_num-1);
cn_num=size(cndata,1);
EndSav=1;
ncp=[0:0.1:1.5]; ncp(12:16)=[1.5 2 3 4 5]; ntw=ncp';
ntw2=[0.5 1 1.5 2 3 4]; ntw2=ntw2';
nmc=1000;
nvmc=zeros(16,nmc,4);
nv=zeros(16,14*4);
sccmc=zeros(15,nmc,4);
scc=zeros(15,14*4);
explist=[3 2 1 4];%tau=5 10 20&30

for exp=1:4
for mc=1:nmc
T=explist(exp);
load(strcat('..\nuclearmonte\output_nuclear_EndSav',num2str(EndSav),'_MC',num2str(mc+(T-1)*1000),'.dat'),'-mat');
dW=zeros(16,2);
for s=1:16
    dW(s,1)=dW(s,1)+sum(output_nuclear(1:(cn_num-1),s,2),1); % global Gt CO2 abated by nuclear power
    for cn=1:(cn_num-1)
        dW(s,2)=dW(s,2)+(output_nuclear(cn,s,3)-output_nuclear(cn,1,3))*output_pop(400,cn); % consumption trillion $
    end
end
nvmc(:,mc,exp)=dW(:,2);
for s=2:15
            A=dW(s-1:s+1,1:2); A=A';
            [rs,ms,bs] = regression(A(1,:),A(2,:));
            sccmc(s,mc,exp)=ms*1000;
end
    A1=dW(1:10,1:2); A1=A1';
    [rs2,ms2,bs2] = regression(A1(1,:),A1(2,:));
    sccmc2(1,mc,exp)=ms2*1000;
    A2=dW(11:16,1:2); A2=A2';
    [rs3,ms3,bs3] = regression(A2(1,:),A2(2,:));
    sccmc2(2,mc,exp)=ms3*1000;    
end
nv(:,(exp-1)*14+1)=mean(nvmc(:,1:nmc,exp),2);
scc(:,(exp-1)*14+1)=mean(sccmc(:,1:nmc,exp),2);
y=prctile(nvmc(:,1:nmc,exp),[1 5 10 20 30 40 50 60 70 80 90 95 99],2);
y2=prctile(sccmc(:,1:nmc,exp),[1 5 10 20 30 40 50 60 70 80 90 95 99],2);
for i=13:-1:1
    nv(:,(exp-1)*14+1+i)=y(:,i);
    scc(:,(exp-1)*14+1+i)=y2(:,i);
end
end
save('..\nuclear\nv.dat','nv');
save('..\nuclear\scc.dat','scc');
