% Script to calc spike stats (fcn of time), for OBvar manuscript
% for a fixed odor 
% calls getSpStats_TimeVary.m

%access to processDataSelectShift.m, getMetaData_OdorDrug_SngleRat_fn.m 
%   exclude_bad_units.m, extractCountsOvrlap.m
% cd ..
% basedir = pwd;
% addpath([basedir '/code_data_analysis/']); 
% cd codeData_CorrTime %get back to current directory

% ODOR: Ethyl Butyrate
odor_to_keep = [2];

prestr = pwd;
% -- consider all 11 rats --
ratst_cell=cell(11,1);
ratst_cell{1}='011416';
ratst_cell{2}='012516'; %used in BGSL
ratst_cell{3}='020416'; %lower responses
ratst_cell{4}='020916';
ratst_cell{5}='021116'; %not as responsive
ratst_cell{6}='021816';
ratst_cell{7}='040115';
ratst_cell{8}='040515';
ratst_cell{9}='040915';
ratst_cell{10}='041615';
ratst_cell{11}='042115';

% Parameters you may want to change
Twin = 0.1; %this has to divide 2 evenly
SponL = 3;
Toverlap = 1; % 0=disjoint, 1=hlfOverlapping
% OUTPUTS--- combine data over ALL 11 rats
mn_allOr = [];
vr_allOr = [];
cov_allOr = [];
crr_eachOr = [];
denom_Or = [];
fano_eachOr = [];

% Compute time-dependent
%     setParam:    Info for calc
%       setParam.Twin:            Time window used for PSTH
%       setParam.SponSampLength: Length of time before stim to use for
%               sampling spontaneous activity
%
setParam = [];
setParam.Twin       = Twin;
setParam.SponSampLength =SponL;
setParam.Toverlap=Toverlap;

nRats=length(ratst_cell);
rats_to_process = 1:nRats;


for ind_Rat=rats_to_process
    %Getting all files in directory
    datfilestr = dir([prestr ratst_cell{ind_Rat} '/' '*.mat']); %Getting all files in directory
    
    % Get metadata
    metaData = getMetaData_OdorDrug_SngleRat_fn(datfilestr,[prestr ratst_cell{ind_Rat}]);
    
    datfilelist = {};
    for j1=1:length(datfilestr)
        datfilelist = [datfilelist [prestr ratst_cell{ind_Rat} '/' datfilestr(j1).name]];
    end
    
    % Now restrict file list to ONLY keep desired files
    HasOdor = sum(metaData(:,2)==odor_to_keep,2);
    HasDrug = sum(metaData(:,3)==0,2);
    
    wantedFiles = find(HasOdor & HasDrug);
    
    datfilelist = datfilelist(wantedFiles);
    
    % Do all the files in the list!
    nFiles = length(datfilelist);
    dat_to_process = 1:nFiles;
    
    for j1=dat_to_process
        [mn_Orth,vr_Orth,cov_Orth,crr_Orth,denmCrr_Orth,~,~,~,~,~,numSpon,numEvok]=...
            getSpStats_TimeVary(datfilelist{j1},setParam);
        
        mn_allOr = [mn_allOr;mn_Orth];
        vr_allOr = [vr_allOr;vr_Orth];
        cov_allOr=[cov_allOr;cov_Orth];
        crr_eachOr=[crr_eachOr;crr_Orth];
        denom_Or=[denom_Or;denmCrr_Orth];
        fanTemp=vr_allOr./mn_allOr;
        fanTemp(isnan(fanTemp))=0; %all with mean=0 set FF to 0
        fano_eachOr=[fano_eachOr; fanTemp];
        
    end

end %for-loop over all rats

tm_spon=Twin*(1:numSpon)'-Twin*numSpon;
tm_evk=(1:numEvok)'*Twin;
tme = [tm_spon;tm_evk];
LnmW=length(tme);

%get FF and Corr pieces via Lin Reg
fano_allOr=zeros(LnmW,1);  fano_allOr_Bic=zeros(LnmW,1); fano_allOr_Mus=zeros(LnmW,1);
crr_allOr=zeros(LnmW,1);  crr_allOr_Bic=zeros(LnmW,1); crr_allOr_Mus=zeros(LnmW,1);
for j=1:LnmW
    fano_allOr(j,1)=mn_allOr(:,j)'*vr_allOr(:,j)/sum(mn_allOr(:,j).^2); %lin regress slope with 0 intercept
    crr_allOr(j,1)=denom_Or(:,j)'*cov_allOr(:,j)/sum(denom_Or(:,j).^2);
end


%% create vars & save pop avgs, store ALL results in expData_OrthEB.mat
psthOr=mean(mn_allOr./Twin);
std_psOr=std(mn_allOr./Twin);

varOr=mean(vr_allOr);
std_vrOr=std(vr_allOr);

FF_Or=mean(fano_eachOr);
std_ffOr=std(fano_eachOr);

CovOr=mean(cov_allOr);
std_cvOr=std(cov_allOr);

corrOr=mean(crr_eachOr);
std_crOr=std(crr_eachOr);

%store all, only -2 to 2s
meanNrnAll=mn_allOr(:,11:end);
varNrnAll=vr_allOr(:,11:end);
covNrnAll=cov_allOr(:,11:end);
denmStd=denom_Or(:,11:end);

save expData_OrthEB tme meanNrnAll varNrnAll covNrnAll denmStd psthOr* varOr* FF_Or* CovOr* corrOr* std_*

disp(['Data is processed and saved in expData_OrthEB.mat. Next parts will plot figures (press any key)'])

pause

% ----- plot results-----
% 2 sub-sections below 

%% plots of Shew data, population averages
load expData_OrthEB.mat 

sclStd=0.2; %scale std deviation, for visual purposes
clGry=.6*ones(1,3);%[1 0 1];
flagEB=1; %1=show error bars, 0=done

figure
subplot(2,1,1)
plot(tme,psthOr,'k','LineWidth',2)
hold on
if(flagEB)
plot(tme,psthOr+sclStd*std_psOr,'-',tme,psthOr-sclStd*std_psOr,'-','color',clGry)
end
set(gca,'FontSize',20)
box off
set(gca,'XLim',[-2+eps 2])
xlabel('Time (s)')
ylabel('Firing Rate (Hz)')

figure
subplot(2,1,1)
plot(tme,varOr,'k','LineWidth',2)
hold on
if(flagEB)
plot(tme,varOr+sclStd*std_vrOr,'-',tme,varOr-sclStd*std_vrOr,'-','color',clGry)
end
set(gca,'FontSize',20)
box off
set(gca,'XLim',[-2+eps 2])
xlabel('Time (s)')
ylabel('Spike Count Variance')

figure
subplot(2,1,1)
plot(tme,FF_Or,'k','LineWidth',2)
hold on
if(flagEB)
plot(tme,FF_Or+sclStd*std_ffOr,'-',tme,FF_Or-sclStd*std_ffOr,'-','color',clGry)
end
set(gca,'FontSize',20)
box off
set(gca,'XLim',[-2+eps 2])
xlabel('Time (s)')
ylabel('Spike Count FF')

figure
subplot(2,1,1)
plot(tme,CovOr,'k','LineWidth',2)
hold on
if(flagEB)
plot(tme,CovOr+sclStd*std_cvOr,'-',tme,CovOr-sclStd*std_cvOr,'-','color',clGry)
end
set(gca,'FontSize',20)
box off
set(gca,'XLim',[-2+eps 2])
xlabel('Time (s)')
ylabel('Spike Count Variance')

figure
subplot(2,1,1)
plot(tme,corrOr,'k','LineWidth',2)
hold on
if(flagEB)
plot(tme,corrOr+sclStd*std_crOr,'-',tme,corrOr-sclStd*std_crOr,'-','color',clGry)
end
set(gca,'FontSize',20)
box off
set(gca,'XLim',[-2+eps 2])
xlabel('Time (s)')
ylabel('Spike Count Correlation')

%% plots of Shew data details, individual cells/pairs
load expData_OrthEB.mat 

iDiv=20; %Twin=.1, 20=2secs
iEv=40;
mnMtch_V=4.5881;%dividing pt
mnMtch_C=7.3627;%dividing pt
   
%!!assuming Twin=100ms!! so 40 total windows!!
msp=mean(meanNrnAll(:,1:iDiv),2);
mev=mean(meanNrnAll(:,iDiv+1:iEv),2);
vsp=mean(varNrnAll(:,1:iDiv),2);
vev=mean(varNrnAll(:,iDiv+1:iEv),2);
csp=mean(covNrnAll(:,1:iDiv),2);
cev=mean(covNrnAll(:,iDiv+1:iEv),2);
denm_sp=mean(denmStd(:,1:iDiv),2);
denm_ev=mean(denmStd(:,iDiv+1:iEv),2);

disp(['Perc of Evoked Var larger than Spont: ',num2str(sum(vev>vsp)/length(vsp))])
disp(['Perc of Evoked Cov larger than Spont: ',num2str(sum(cev>csp)/length(csp))])

idz=vsp~=0;
vsp=vsp(idz);
vev=vev(idz);
msp=msp(idz);
mev=mev(idz);
figure
for j=1:length(msp)
    semilogy([msp(j) mev(j)],'.-','color',.5*ones(1,3),'MarkerSize',12)
    if(j==1)
        hold on
    end
end
set(gca,'FontSize',24)
axis([.8 2.2 0 14])
ylabel('Mean')
semilogy(mean([msp mev]),'ks-','LineWidth',2,'MarkerSize',12,'MarkerFaceColor','k')

figure
for j=1:length(vsp)
    semilogy([vsp(j) vev(j)],'.-','color',.5*ones(1,3),'MarkerSize',12)
    if(j==1)
        hold on
    end
end
set(gca,'FontSize',24)
axis([.8 2.2 .0017 51.1416])
ylabel('Var')
%pause
semilogy(mean([vsp vev]),'ks-','LineWidth',2,'MarkerSize',12,'MarkerFaceColor','k')

idz=csp~=0;
csp=csp(idz);
cev=cev(idz);
denm_sp=denm_sp(idz);
denm_ev=denm_ev(idz);
figure
hold on
for j=1:length(csp)
    plot([csp(j) cev(j)],'.-','color',.5*ones(1,3),'MarkerSize',12)
end
set(gca,'FontSize',24)
axis([.8 2.2 -1 3.5])
ylabel('Covar')
%pause
plot(mean([csp cev]),'ks-','LineWidth',2,'MarkerSize',12,'MarkerFaceColor','k')  

FFsp_reg=msp'*vsp/sum(msp.^2); %FanoFact from linear regress
idSp=msp<mnMtch_V;
FFsp2=msp(idSp)'*vsp(idSp)/sum(msp(idSp).^2);
FFev_reg=mev'*vev/sum(mev.^2); %FanoFact from linear regress
idEv=mev<mnMtch_V;
FFev2=mev(idEv)'*vev(idEv)/sum(mev(idEv).^2);
figure
hold on
box off
plot(msp,vsp,'ko','MarkerFaceColor','k')
    plot(msp,msp*FFsp_reg,'k','LineWidth',1)
    plot(msp(idSp),msp(idSp)*FFsp2,'color',.5*ones(1,3),'LineWidth',1)
plot(mev,vev,'ro')
    plot(mev,mev*FFev_reg,'r','LineWidth',1)
    plot(mev,mev*FFev2,'r','LineWidth',1)
set(gca,'FontSize',24)
xlabel('Mean Spikes')
ylabel('Variance')
axis([0 13.5 0 35])
% -- re-do linear regres to get stats --
%[fit1,gof,fitinfo] = fit(msp,vsp,'poly1');
%lmVsp=fitlm(msp,vsp); lmVsp
%[ms,~,]=lsqlin([ones(length(msp),1) msp],vsp,[],[],[],[],zeros(2,1),[0 10]);
Rsq_Vsp=sum((FFsp_reg*msp-mean(vsp)).^2)/sum((vsp-mean(vsp)).^2)
Rsq_Vev=sum((FFev_reg*mev-mean(vev)).^2)/sum((vev-mean(vev)).^2)
Rsq_Vev2=sum((FFev2*mev(idEv)-mean(vev(idEv))).^2)/sum((vev(idEv)-mean(vev(idEv))).^2)

corrSp_reg=denm_sp'*csp/sum(denm_sp.^2);
    idSp=denm_sp<mnMtch_C;
corrSp2=denm_sp(idSp)'*csp(idSp)/sum(denm_sp(idSp).^2);
corrEv_reg=denm_ev'*cev/sum(denm_ev.^2);
    idEv=denm_ev<mnMtch_C;
corrEv2=denm_ev(idEv)'*cev(idEv)/sum(denm_ev(idEv).^2);
figure
hold on
box off
plot(denm_sp,csp,'ko','MarkerFaceColor','k')
    plot(denm_sp,denm_sp*corrSp_reg,'k','LineWidth',1)
    plot(denm_sp(idSp),denm_sp(idSp)*corrSp2,'color',.5*ones(1,3),'LineWidth',1)
plot(denm_ev,cev,'ro')
    plot(denm_ev,denm_ev*corrEv_reg,'r','LineWidth',1)
    plot(denm_ev,denm_ev*corrEv2,'r','LineWidth',1)
set(gca,'FontSize',24)
xlabel('Sqrt')
ylabel('Covariance')   
% get R^2 
Rsq_Csp=1-sum((corrSp_reg*denm_sp-csp).^2)/sum((csp-mean(csp)).^2)
Rsq_Cev=1-sum((corrEv_reg*denm_ev-cev).^2)/sum((cev-mean(cev)).^2)
Rsq_Cev2=1-sum((corrEv2*denm_ev(idEv)-cev(idEv)).^2)/sum((cev(idEv)-mean(cev(idEv))).^2)
