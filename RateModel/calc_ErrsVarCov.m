%script to results for AWAKE data (Bolding,Franks18)
% AFTER getting best 100, calc errors w.r.t Var & Cov HERE 
% (run calc_ErrsPSTH.m first to get dCovAwk2.mat)
% see calc_AwakeErrs.m for PSTH, 
% assuming have sSm[j].mat, sVry[j].mat, sLrg[j].mat for j=1,..,100 FROM mcRateAwk.m
% assuming have frWorse2.mat that has WORSE firing rate model results (gray curves in figure)

cc=[zeros(1,3); [189 130 0]/255 ; 0 0 1];       %solid color
cc2=[51 51 51; 242 230 204 ; 204 204 255]./255; %transparent color
    
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

Twin=0.1;

numGs=100;

% --- outputs --- 
l1w_Var=zeros(numGs,3); %L-1 errors, NORMALIZD; 1st col is _sc, 2nd col is _n (sigVary) 
nl1_Var=zeros(numGs,3); %errors WITH signs
indWrs_Var=zeros(2,3); %2 worse for Var

l1w_Cov=zeros(numGs,3); %L-1 errors, NORMALIZD; 1st col is _sc, 2nd col is _n (sigVary) 
nl1_Cov=zeros(numGs,3); %errors WITH signs
indWrs_Cov=zeros(2,3); %2 worse for cov

l1e_Fr=zeros(numGs,3); %L-1 errors, NORMALIZD

% for geting L1-error
load('expD_awake_EB.mat','tme','varOr','CovOr','psthOr') %experim data
dataPs=interp1(tme',psthOr,tmv','pchip');
iSt_D=round((-0.4-tme(1))/Twin)+1;
iEn_D=round((0.5-tme(1))/Twin)+1;
%get L1-error
for j=1:numGs
    load(['sSm',num2str(j)],'varFRt_tw','covFRt_tw','mcFm')
    l1w_Var(j,1)=sum(abs(varOr(iSt_D:iEn_D)-varFRt_tw))/sqrt(sum(varOr(iSt_D:iEn_D).^2));
    nl1_Var(j,1)=sum(varOr(iSt_D:iEn_D)-varFRt_tw)*Twin;
    l1w_Cov(j,1)=sum(abs(CovOr(iSt_D:iEn_D)-covFRt_tw))/sqrt(sum(CovOr(iSt_D:iEn_D).^2));
    nl1_Cov(j,1)=sum(CovOr(iSt_D:iEn_D)-covFRt_tw)*Twin;
    l1e_Fr(j,1)=sum(abs(dataPs-mcFm))/sqrt(sum(dataPs.^2));
    
    load(['sVry',num2str(j)],'varFRt_tw','covFRt_tw','mcFm')
    l1w_Var(j,2)=sum(abs(varOr(iSt_D:iEn_D)-varFRt_tw))/sqrt(sum(varOr(iSt_D:iEn_D).^2));
    nl1_Var(j,2)=sum(varOr(iSt_D:iEn_D)-varFRt_tw)*Twin;
    l1w_Cov(j,2)=sum(abs(CovOr(iSt_D:iEn_D)-covFRt_tw))/sqrt(sum(CovOr(iSt_D:iEn_D).^2));
    nl1_Cov(j,2)=sum(CovOr(iSt_D:iEn_D)-covFRt_tw)*Twin;
    l1e_Fr(j,2)=sum(abs(dataPs-mcFm))/sqrt(sum(dataPs.^2));
    
    load(['sLrg',num2str(j)],'varFRt_tw','covFRt_tw','mcFm')
    l1w_Var(j,3)=sum(abs(varOr(iSt_D:iEn_D)-varFRt_tw))/sqrt(sum(varOr(iSt_D:iEn_D).^2));
    nl1_Var(j,3)=sum(varOr(iSt_D:iEn_D)-varFRt_tw)*Twin;
    l1w_Cov(j,3)=sum(abs(CovOr(iSt_D:iEn_D)-covFRt_tw))/sqrt(sum(CovOr(iSt_D:iEn_D).^2));
    nl1_Cov(j,3)=sum(CovOr(iSt_D:iEn_D)-covFRt_tw)*Twin;
    l1e_Fr(j,3)=sum(abs(dataPs-mcFm))/sqrt(sum(dataPs.^2));
end

[m,i]=min(nl1_Var);
indWrs_Var(1,:)=i;
[m,i]=max(nl1_Var);
indWrs_Var(2,:)=i;

[m,i]=min(nl1_Cov);
indWrs_Cov(1,:)=i;
[m,i]=max(nl1_Cov);
indWrs_Cov(2,:)=i;


% determine tenBst here:
tenBst=zeros(10,3); 
TotErr=l1e_Fr+l1w_Var+l1w_Cov;
for j=1:3
    [m,i]=sort(TotErr(:,j));
    tenBst(:,j)=i(1:10);
end

%for plotting 10 best
tM=(-.4:.1:.5)';

figure %PSTH here
for k=1:3
    %[m,ind]=sort(l1e_Fr(:,k));
    ind=tenBst(:,k);
    subplot(1,3,k)
    switch k
        case 1
            flN='sSm';
        case 2
            flN='sVry';
        case 3
            flN='sLrg';
    end
    hold on
    for j=1:10
        load([flN,num2str(ind(j))])
        plot(tmv',mcFm,'b')

    end
    load frWorse2.mat
    plot(tmv,mcFrt_Worse(2*k-1,:),'--','color',.4*ones(1,3))
    plot(tmv,mcFrt_Worse(2*k,:),'--','color',.4*ones(1,3))
    plot(tmv,dataPs,'k','LineWidth',3)
    set(gca,'FontSize',18)
    set(gca,'XLim',[-.5 .5])
end

figure %do var here
for k=1:3
    %[m,ind]=sort(l1w_Var(:,k));
    ind=tenBst(:,k);
    subplot(1,3,k)
    switch k
        case 1
            flN='sSm';
        case 2
            flN='sVry';
        case 3
            flN='sLrg';
    end
    hold on
    for j=1:10
        load([flN,num2str(ind(j))])
        plot(tM,varFRt_tw,'b')

    end
    load([flN,num2str(indWrs_Var(1,k))])
    plot(tM,varFRt_tw,'--','color',.4*ones(1,3))
    load([flN,num2str(indWrs_Var(2,k))])
    plot(tM,varFRt_tw,'--','color',.4*ones(1,3))
    plot(tme,varOr,'k','LineWidth',3)
    set(gca,'FontSize',18)
    set(gca,'XLim',[-.5 .5])
end

figure %now do Cov
for k=1:3 
    ind=tenBst(:,k);
    subplot(1,3,k)
    switch k
        case 1
            flN='sSm';
        case 2
            flN='sVry';
        case 3
            flN='sLrg';
    end
    hold on
    for j=1:10
        load([flN,num2str(ind(j))])
        plot(tM,covFRt_tw,'b')

    end
    load([flN,num2str(indWrs_Cov(1,k))])
    plot(tM,covFRt_tw,'--','color',.4*ones(1,3))
    load([flN,num2str(indWrs_Cov(2,k))])
    plot(tme,CovOr,'k','LineWidth',3)
    set(gca,'FontSize',18)
    set(gca,'XLim',[-.5 .5])
end

%% show 10 best coupling strengths, assuming indexes are in tenBst(10,3)
load dCovAwk2.mat
figure('Renderer', 'Painters');
hold on
plot3(gG2Mv(indCT_sc(tenBst(:,1))),gM2Gv(indCT_sc(tenBst(:,1))),gGCIv(indCT_sc(tenBst(:,1))),'.','color',cc2(1,:),'MarkerFaceColor',cc2(1,:),'MarkerSize',50)

plot3(gG2Mv(indCT_n(tenBst(:,2))),gM2Gv(indCT_n(tenBst(:,2))),gGCIv(indCT_n(tenBst(:,2))),'.','color',cc2(2,:),'MarkerFaceColor',cc2(2,:),'MarkerSize',50)

plot3(gG2Mv(indCT_cs(tenBst(:,3))),gM2Gv(indCT_cs(tenBst(:,3))),gGCIv(indCT_cs(tenBst(:,3))),'.','color',cc2(3,:),'MarkerFaceColor',cc2(3,:),'MarkerSize',50)
grid on
view([-43 34])
xlabel('G to M'); ylabel('M to G');
axis([0 10 0 10 0 5])

% displays avg coupling strength for Table 3
A_sc=[gG2Mv(indCT_sc(tenBst(:,1)))  gGCIv(indCT_sc(tenBst(:,1))) gM2Gv(indCT_sc(tenBst(:,1)))];
mean(A_sc)
A_n=[gG2Mv(indCT_n(tenBst(:,2))) gGCIv(indCT_n(tenBst(:,2))) gM2Gv(indCT_n(tenBst(:,2))) ];
mean(A_n)
A_cs=[gG2Mv(indCT_cs(tenBst(:,3))) gGCIv(indCT_cs(tenBst(:,3))) gM2Gv(indCT_cs(tenBst(:,3))) ];
mean(A_cs)


%% box plots for errors for 10 best, assuming indexes are in ?tenBst(10,3)?
sErr_C_sc=l1w_Cov(tenBst(:,1),1);
sErr_C_n=l1w_Cov(tenBst(:,2),2);
sErr_C_cs=l1w_Cov(tenBst(:,3),3);
sErr_V_sc=l1w_Var(tenBst(:,1),1);
sErr_V_n=l1w_Var(tenBst(:,2),2);
sErr_V_cs=l1w_Var(tenBst(:,3),3);
Scl_TErr_sc=l1e_Fr(tenBst(:,1),1)+sErr_V_sc+sErr_C_sc;
Scl_TErr_n=l1e_Fr(tenBst(:,2),2)+sErr_V_n+sErr_C_n;
Scl_TErr_cs=l1e_Fr(tenBst(:,3),3)+sErr_V_cs+sErr_C_cs;
figure
hold on
boxchart(ones(length(Scl_TErr_sc),1),Scl_TErr_sc,'MarkerStyle','none','BoxFaceColor',cc(1,:),'WhiskerLineColor',cc(1,:))
boxchart(2*ones(length(Scl_TErr_n),1),Scl_TErr_n,'MarkerStyle','none','BoxFaceColor',cc(2,:),'WhiskerLineColor',cc(2,:))
boxchart(3*ones(length(Scl_TErr_cs),1),Scl_TErr_cs,'MarkerStyle','none','BoxFaceColor',cc(3,:),'WhiskerLineColor',cc(3,:))
set(gca,'FontSize',18)
ylabel('Total Error')

% var
figure
hold on
boxchart(ones(length(Scl_TErr_sc),1),sErr_V_sc,'MarkerStyle','none','BoxFaceColor',cc(1,:),'WhiskerLineColor',cc(1,:))
boxchart(2*ones(length(Scl_TErr_n),1),sErr_V_n,'MarkerStyle','none','BoxFaceColor',cc(2,:),'WhiskerLineColor',cc(2,:))
boxchart(3*ones(length(Scl_TErr_cs),1),sErr_V_cs,'MarkerStyle','none','BoxFaceColor',cc(3,:),'WhiskerLineColor',cc(3,:))
set(gca,'FontSize',18)
ylabel('Error of Var')

% Covar
figure
hold on
boxchart(ones(length(Scl_TErr_sc),1),sErr_C_sc,'MarkerStyle','none','BoxFaceColor',cc(1,:),'WhiskerLineColor',cc(1,:))
boxchart(2*ones(length(Scl_TErr_n),1),sErr_C_n,'MarkerStyle','none','BoxFaceColor',cc(2,:),'WhiskerLineColor',cc(2,:))
boxchart(3*ones(length(Scl_TErr_cs),1),sErr_C_cs,'MarkerStyle','none','BoxFaceColor',cc(3,:),'WhiskerLineColor',cc(3,:))
set(gca,'FontSize',18)
ylabel('Error of Covar')

