% % Author: Rong Wang // contact rongwang@fudan.edu.cn //
% % Date: 2022.9.24
% % 
% 
hist=[279 160 142 29 26 8 1 2];
N=sum(hist(2:8),2);
randnu=zeros(1000,6+10*N+134+1000);
M=randperm(N);
W=zeros(N,1);
k=0;
for i=1:7
    for j=1:hist(i+1)
        k=k+1;
        W(M(k),1)=i;
    end
end
for i=1:1000
    for j=1:5
        randnu(i,j)=0.33+0.33*rand; % threshold of probability for tipping 1-5
    end
    randnu(i,6)=min(2,max(-2,randn)); % threshold of cumulative nuclear PWh for accident
    for j=1:1000
        randnu(i,6+10*N+134+j)=randn;% threshold of cumulative nuclear PWh for accident
        while randnu(i,6+10*N+134+j)>2 || randnu(i,6+10*N+134+j)<-2
            randnu(i,6+10*N+134+j)=randn;
        end
    end

    for j=1:(N*10)
        randnu(i,j+6)=W(mod(j-1,N)+1,1);
    end
    N2=randperm(N*10);
    randnu(i,(N*10+7):(N*10+7+133))=N2(1:134); % start point: uniform distribution 1-7
end

save('..\nuclearmonte\randnu.dat','randnu');



% % for ECS
% esc2=zeros(10000,1);
% for i=1:10000
%     esc2(i) = 3.25 + 1.3/4 * randn; % (3.9-2.6)/2/2=0.32 as sigma
% %     esc2(i) = min(3.9,max(2.6,3.25 + 1.3/4 * randn)); % (3.9-2.6)/2/2=0.32 as sigma
% end
% yy=prctile(esc2,[10:10:90]);
% yy(10)=2.6; yy(11)=3.9; yy(12)=3.1;

