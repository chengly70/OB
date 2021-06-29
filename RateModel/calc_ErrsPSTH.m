%script to get PSTH results for AWAKE data (Bolding,Franks18)
% for anesthetized, very similar 
% saves indices for 100 best (PSTH) in dCovAwk2.mat
% saves indices for 2 worst (PSTH) in dAwk_Worse2.mat

cc=[zeros(1,3); [189 130 0]/255 ; 0 0 1];
    
load dSeqHalton2 %gets nvcSmp [0,1] in 3 dims
gP2M=-5;
gM2P=5;
gG2Mv=-5*2*nvcSmp(:,1);
gGCIv=-2.5*2*nvcSmp(:,2);
gM2Gv=5*2*nvcSmp(:,3);
%make it positive just for ease
gGCIv=abs(gGCIv);
gG2Mv=abs(gG2Mv);
    
dt=0.001;
tmv=(-0.5 : dt : 0.5)'; %in seconds
Lt=length(tmv); 

numGs=10000;

% --- outputs --- 
l1w_err=zeros(numGs,3); %L-1 errors; 1st col is _sc, 2nd col is _n (sigVary)
nl1_err=zeros(numGs,3); %errors WITH signs
indWorse=zeros(2,3); %2 worse PSTH curves for each

% for geting L1-error
load('expD_awake_EB.mat','tme','psthOr','std_psOr') %experim data
dataPs=interp1(tme',psthOr,tmv','pchip');
datStdPs=interp1(tme',std_psOr,tmv','pchip');

%get calc data (MC in fr model)
load(['frt2_sc.mat'],'mcFrate')
mcFr_sc=mcFrate;
load(['frt2_n.mat'],'mcFrate')
mcFr_n=mcFrate;
load(['frt2_cs.mat'],'mcFrate')
mcFr_cs=mcFrate;

%get L1-error
iStrt=round((-.2-tmv(1))/dt)+1;
iEnd=round((.2-tmv(1))/dt)+1;
for j=1:numGs
    l1w_err(j,1)=sum(abs(dataPs(iStrt:iEnd)-mcFr_sc(j,iStrt:iEnd)))*dt;
    l1w_err(j,2)=sum(abs(dataPs(iStrt:iEnd)-mcFr_n(j,iStrt:iEnd)))*dt;
    l1w_err(j,3)=sum(abs(dataPs(iStrt:iEnd)-mcFr_cs(j,iStrt:iEnd)))*dt;
    
    nl1_err(j,1)=sum((dataPs(iStrt:iEnd)-mcFr_sc(j,iStrt:iEnd)))*dt;
    nl1_err(j,2)=sum((dataPs(iStrt:iEnd)-mcFr_n(j,iStrt:iEnd)))*dt;
    nl1_err(j,3)=sum((dataPs(iStrt:iEnd)-mcFr_cs(j,iStrt:iEnd)))*dt;
end

[m,i]=min(nl1_err);
indWorse(1,:)=i;
[m,i]=max(nl1_err);
indWorse(2,:)=i;

ptSize=40;

% get absolute best fits
[m,in1]=min(l1w_err(:,1)); [m,in2]=min(l1w_err(:,2)); [m,in3]=min(l1w_err(:,3));

% display coupling strengths for best fit, not used in paper just FYI
figure
hold on
plot(tmv,mcFr_sc(in1,:),'color',cc(1,:))
plot(tmv,mcFr_n(in2,:),'color',cc(2,:))
plot(tmv,mcFr_cs(in3,:),'color',cc(3,:))
plot(tme,psthOr,'b','LineWidth',3)
set(gca,'FontSize',18)
set(gca,'XLim',[-0.5 0.5])
[gG2Mv(in1), gGCIv(in1), gM2Gv(in1)]
[gG2Mv(in2), gGCIv(in2), gM2Gv(in2)]
[gG2Mv(in3), gGCIv(in3), gM2Gv(in3)]


%suppress indiv cell results
    
    for j=1:3
        figure
        hold on
        %scatter3(gG2Mv,gM2Gv,gGCIv,ptSize,col2,'fill');
        gG2MvA=gG2Mv(l1w_err(:,j)<thresErr); gGCIvA=gGCIv(l1w_err(:,j)<thresErr); gM2GvA=gM2Gv(l1w_err(:,j)<thresErr); colA=zeros(length(gM2GvA),1);
        gG2MvB=gG2Mv(l1w_err(:,j)>=thresErr); gGCIvB=gGCIv(l1w_err(:,j)>=thresErr); gM2GvB=gM2Gv(l1w_err(:,j)>=thresErr); colB=ones(length(gM2GvB),1);
        scatter3(gG2MvB,gM2GvB,gGCIvB,ptSize,colB,'o');
        scatter3(gG2MvA,gM2GvA,gGCIvA,ptSize,colA,'filled');
        set(gca,'FontSize',18)
        xlabel('G2M')
        ylabel('M2G')
        zlabel('G common')
        grid on
        colormap jet
        axis([0 10 0 10 0 5])
        view([-43 34])
    end
% %% save results in dCovAwk.mat; WITHIN Threshold
% ind=(1:10000)';
% indCT_sc=ind(colTot(:,1)==0); indCT_n=ind(colTot(:,2)==0); indCT_cs=ind(colTot(:,3)==0);
% save dCovAwk indCT_n indCT_sc indCT_cs

%% save results in dCovAwk.mat, just save top 100
ind=(1:10000)';
[m,i]=sort(l1w_err(:,1)); indCT_sc=i(1:100);
[m,i]=sort(l1w_err(:,2)); indCT_n=i(1:100);
[m,i]=sort(l1w_err(:,3)); indCT_cs=i(1:100);
save dCovAwk2 indCT_n indCT_sc indCT_cs

%%

figure
for k=1:3
    
    switch k
        case 1
            mcFr=mcFr_sc;
            indThr=indCT_sc;
        case 2
            mcFr=mcFr_n;
            indThr=indCT_n;
        case 3
            mcFr=mcFr_cs;
            indThr=indCT_cs;
    end
    subplot(1,3,k)
    hold on
    for j=1:10
        
        plot(tmv,mcFr(indThr(j),:),'b')
        
    end
    plot(tmv,mcFr(indWorse(1,k),:),'--','color',.4*ones(1,3))
    plot(tmv,mcFr(indWorse(2,k),:),'--','color',.4*ones(1,3))
    plot(tme,psthOr,'k','LineWidth',3)
    set(gca,'FontSize',18)
    set(gca,'XLim',[-.5 .5])
end

% save the worse ones in dAwk_Worse2
save dAwk_Worse2 indWorse

% shows avg coupling strength for best 100 (not shown in paper)
load dCovAwk2
mean([gG2Mv(indCT_sc), gM2Gv(indCT_sc), gGCIv(indCT_sc)])
mean([gG2Mv(indCT_n), gM2Gv(indCT_n), gGCIv(indCT_n)])
mean([gG2Mv(indCT_cs), gM2Gv(indCT_cs), gGCIv(indCT_cs)])
