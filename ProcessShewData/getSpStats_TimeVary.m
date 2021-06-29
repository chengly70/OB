function [mn_Orth,vr_Orth,cov_Orth,crr_Orth,denmCrr_Orth,mn_Retr,vr_Retr,cov_Retr,crr_Retr,denmCrr_Retr,numSpon,numEvok,mn_spkCn,vr_spkCn,cov_spkCn,crr_spkCn,denmCrr]...
    =getSpStats_TimeVary(whichDatS,setParam)
%% getSpStats_TimeVary: function to produce stats
% 
% INPUTS
%   whichDatS:    ONLY support file name
%   setParam:     Info for selectivity calc
%       setParam.Twin:  Time window used for PSTH
%       setParam.SponSampLength: Length of time before stim to use for
%               sampling spontaneous activity

% OUTPUTS
%   [mn_Orth,vr_Orth,cov_Orth,crr_Orth,denmCrr_Orth,mn_Retr,vr_Retr,cov_Retr,crr_Retr,denmCrr_Retr,numSpon,numEvok,mn_spkCn,vr_spkCn,cov_spkCn,crr_spkCn,denmCrr]
% whole times series (spon&Evok) of mean/var/cov/corr spike counts in
% Orth/Retr, respectively.  numSpon & numEvok determine state in whole time series
%


% Time window for samples
Twin   = setParam.Twin;

%length of time to sample before stim
SponSampLength = setParam.SponSampLength;

WinOvrLap=setParam.Toverlap; %0=disjoint windows, 1=overlapping windows
% numEvok is an output
numEvok = (2/Twin); %assuming non-overlapping windows

%% Shift each trial by 15 seconds for display and calculation
StimShift = 15;
FirstEvok = round(StimShift/Twin) + 1;
LastEvok  = FirstEvok+numEvok-1;

% run processDataSelectShift.m
[stim_trials,nEpoch,nOB,~]=processDataSelectShift(Twin,0,StimShift,whichDatS); %always set WinOvrLap=0 in processDataSelectShift b/c extractCountsOvrLap.m is wrong

lenTrial = size(stim_trials{3}.OBcounts,1); %30 sec per trial but some extra time before start of next

OBcount_allTrial_allCell = zeros(lenTrial,nEpoch,nOB);
for k=1:nOB
    OBcount_allTrial = zeros(lenTrial,nEpoch);
    how_many_stim_trials = 0;
    for k1=1:nEpoch
        if (stim_trials{k1}.stim_at_startT==1)
           % Is there a stimulus during this epoch??
           how_many_stim_trials = how_many_stim_trials+1;
        
            % Read data, put into array
            [nwin,nU]=size(stim_trials{k1}.OBcounts);

            if (k<=nU)   % Otherwise, No spikes belonging to this unit on that trial
                         % Recorded array was set based on max UID that
                         % appeared in that trial
                if(nwin<lenTrial)
                    OBcount_allTrial(1:nwin,how_many_stim_trials)=stim_trials{k1}.OBcounts(1:nwin,k);
                else
                    OBcount_allTrial(1:lenTrial,how_many_stim_trials)=stim_trials{k1}.OBcounts(1:lenTrial,k);
                end
            end
        end
    end
    OBcount_allTrial_allCell(:,:,k) = OBcount_allTrial;
end
% Set properly to # of trials with stim
OBcount_allTrial_allCell = OBcount_allTrial_allCell(:,1:how_many_stim_trials,:);
% Use same criteria as before to exclude cells with no spikes 
% Average across trials, for each OB cell
OBmeankeep = [];
for k=1:nOB
    OBcount_allTrial = OBcount_allTrial_allCell(:,:,k);

    OBmeankeep = [OBmeankeep; mean(OBcount_allTrial')];
end
% Exclude firing rates that are identically 0
i_OBs=(sum(OBmeankeep')~=0);    %also used in loop below for Cov/corr
OBmeankeep=OBmeankeep(i_OBs',:);
nOB=sum(i_OBs); %update size of OB
if(nOB ~= size(OBcount_allTrial_allCell,3))
    %throw out some cell(s)
    OBcount_allTrial_allCell=OBcount_allTrial_allCell(:,:,i_OBs);
end

if(WinOvrLap==0) %disjoint windows
    %% Data calculation
    kDat = 1:nOB;
    
    % To take samples of spontaneous activity
    LastSpon  = FirstEvok-1;
    FirstSpon = round((StimShift-SponSampLength)/Twin)+1 ;   %SponSampLength (5) seconds before
    nmBns = (LastEvok-FirstSpon+1); %total number of time bins
    %--- OUTPUTS ---
    numSpon=LastSpon-FirstSpon+1;
    mn_Orth=zeros(nOB,nmBns);
    vr_Orth=zeros(nOB,nmBns);
    cov_Orth=zeros(nOB*(nOB-1)/2,nmBns);
    crr_Orth=zeros(nOB*(nOB-1)/2,nmBns); %padd any with 0 var with 0 corr
    denmCrr_Orth=zeros(nOB*(nOB-1)/2,nmBns);
    mn_Retr=zeros(nOB,nmBns);
    vr_Retr=zeros(nOB,nmBns);
    cov_Retr=zeros(nOB*(nOB-1)/2,nmBns);
    crr_Retr=zeros(nOB*(nOB-1)/2,nmBns); %padd any with 0 var with 0 corr
    denmCrr_Retr=zeros(nOB*(nOB-1)/2,nmBns);
    mn_spkCn=zeros(nOB,nmBns);
    vr_spkCn=zeros(nOB,nmBns);
    cov_spkCn=zeros(nOB*(nOB-1)/2,nmBns);
    crr_spkCn=zeros(nOB*(nOB-1)/2,nmBns); %padd any with 0 var with 0 corr
    denmCrr=zeros(nOB*(nOB-1)/2,nmBns);
    
    nTrial   = size(OBcount_allTrial_allCell,2); %all are 20
    
    
    for k=kDat
        % For each cell
        OBtemp = OBcount_allTrial_allCell(FirstSpon:LastEvok,:,k);
        
        mn_Orth(k,:)=mean(OBtemp(:,1:10)');
        vr_Orth(k,:)=var(OBtemp(:,1:10)');
        mn_Retr(k,:)=mean(OBtemp(:,11:20)');
        vr_Retr(k,:)=var(OBtemp(:,11:20)');
        mn_spkCn(k,:)=mean(OBtemp');
        vr_spkCn(k,:)=var(OBtemp');
        
    end
    
    %indicies for Cov (non-var)
    linInd=[];
    for rwI=1:(nOB-1)
        for clI=rwI+1:nOB
            linInd=[linInd; sub2ind([nOB nOB],rwI,clI)];
        end
    end
    
    %loop over time bins, not too inefficient
    for j=1:nmBns
        % 20 x nOB matrix
        tmpOB=squeeze(OBcount_allTrial_allCell(FirstSpon+(j-1),:,:)); %get all data
        tmpOB=tmpOB-repmat(mn_spkCn(:,j)',nTrial,1); %centered
        covM_tmp=tmpOB'*tmpOB./(nTrial-1); %unbiased estim Cov; nOB x nOB matrix
        tmpVr_trun=vr_spkCn(:,j);
        tmpVr_trun(tmpVr_trun<1e-4)=1e-4; %lower Bound on Var for correl calc
        corM_tmp=covM_tmp./sqrt(tmpVr_trun*tmpVr_trun');
        cov_spkCn(:,j)=covM_tmp(linInd);
        crr_spkCn(:,j)=corM_tmp(linInd);
        corM_tmp=sqrt(tmpVr_trun*tmpVr_trun');
        denmCrr(:,j)=corM_tmp(linInd);
        
        % Ortho: 10 x nOB matrix
        tmpOB=squeeze(OBcount_allTrial_allCell(FirstSpon+(j-1),1:10,:)); %get all data
        tmpOB=tmpOB-repmat(mn_Orth(:,j)',nTrial/2,1); %centered
        covM_tmp=tmpOB'*tmpOB./(nTrial/2-1); %unbiased estim Cov; nOB x nOB matrix
        tmpVr_trun=vr_Orth(:,j);
        tmpVr_trun(tmpVr_trun<1e-4)=1e-4; %lower Bound on Var for correl calc
        corM_tmp=covM_tmp./sqrt(tmpVr_trun*tmpVr_trun');
        cov_Orth(:,j)=covM_tmp(linInd);
        crr_Orth(:,j)=corM_tmp(linInd);
        corM_tmp=sqrt(tmpVr_trun*tmpVr_trun');
        denmCrr_Orth(:,j)=corM_tmp(linInd);
        
        % Retro: 10 x nOB matrix
        tmpOB=squeeze(OBcount_allTrial_allCell(FirstSpon+(j-1),11:20,:)); %get all data
        tmpOB=tmpOB-repmat(mn_Retr(:,j)',nTrial/2,1); %centered
        covM_tmp=tmpOB'*tmpOB./(nTrial/2-1); %unbiased estim Cov; nOB x nOB matrix
        tmpVr_trun=vr_Retr(:,j);
        tmpVr_trun(tmpVr_trun<1e-4)=1e-4; %lower Bound on Var for correl calc
        corM_tmp=covM_tmp./sqrt(tmpVr_trun*tmpVr_trun');
        cov_Retr(:,j)=covM_tmp(linInd);
        crr_Retr(:,j)=corM_tmp(linInd);
        corM_tmp=sqrt(tmpVr_trun*tmpVr_trun');
        denmCrr_Retr(:,j)=corM_tmp(linInd);
    end


elseif(WinOvrLap==1)
    %% do it again on half-grid IF WinOvrLap==1 run processDataSelectShift.m
    [stim_trials,nEpoch,nOB,~]=processDataSelectShift(Twin,0,StimShift+.5*Twin,whichDatS); %always set WinOvrLap=0 in processDataSelectShift b/c extractCountsOvrLap.m is wrong
    
    lenTrial = size(stim_trials{3}.OBcounts,1); %30 sec per trial but some extra time before start of next
    
    OBcount_allTrial_allCellHlf = zeros(lenTrial,nEpoch,nOB);
    for k=1:nOB
        OBcount_allTrial = zeros(lenTrial,nEpoch);
        how_many_stim_trials = 0;
        for k1=1:nEpoch
            if (stim_trials{k1}.stim_at_startT==1)
                % Is there a stimulus during this epoch??
                how_many_stim_trials = how_many_stim_trials+1;
                
                % Read data, put into array
                [nwin,nU]=size(stim_trials{k1}.OBcounts);
                
                if (k<=nU)   % Otherwise, No spikes belonging to this unit on that trial
                    % Recorded array was set based on max UID that
                    % appeared in that trial
                    if(nwin<lenTrial)
                        OBcount_allTrial(1:nwin,how_many_stim_trials)=stim_trials{k1}.OBcounts(1:nwin,k);
                    else
                        OBcount_allTrial(1:lenTrial,how_many_stim_trials)=stim_trials{k1}.OBcounts(1:lenTrial,k);
                    end
                end
            end
        end
        OBcount_allTrial_allCellHlf(:,:,k) = OBcount_allTrial;
    end
    % Set properly to # of trials with stim
    OBcount_allTrial_allCellHlf = OBcount_allTrial_allCellHlf(:,1:how_many_stim_trials,:);
    % Use same criteria as before to exclude cells with no spikes
    % Average across trials, for each OB cell
    OBmeankeep = [];
    for k=1:nOB
        OBcount_allTrial = OBcount_allTrial_allCellHlf(:,:,k);
        
        OBmeankeep = [OBmeankeep; mean(OBcount_allTrial')];
    end
    % Exclude firing rates that are identically 0
    i_OBs=(sum(OBmeankeep')~=0);    %also used in loop below for Cov/corr
    OBmeankeep=OBmeankeep(i_OBs',:);
    nOB=sum(i_OBs); %update size of OB
    if(nOB ~= size(OBcount_allTrial_allCellHlf,3))
        %throw out some cell(s)
        OBcount_allTrial_allCellHlf=OBcount_allTrial_allCellHlf(:,:,i_OBs);
    end


    %% Data calculation
    kDat = 1:nOB;
    
    % To take samples of spontaneous activity
    LastSpon  = FirstEvok-1;
    FirstSpon = round((StimShift-SponSampLength)/Twin)+1 ;   %SponSampLength (5) seconds before
    nmBns = (LastEvok-FirstSpon+1); %total number of time bins
    %--- OUTPUTS ---
    numSpon=LastSpon-FirstSpon+1;
    mn_Orth=zeros(nOB,nmBns);
    vr_Orth=zeros(nOB,nmBns);
    cov_Orth=zeros(nOB*(nOB-1)/2,nmBns);
    crr_Orth=zeros(nOB*(nOB-1)/2,nmBns); %padd any with 0 var with 0 corr
    denmCrr_Orth=zeros(nOB*(nOB-1)/2,nmBns);
    mn_Retr=zeros(nOB,nmBns);
    vr_Retr=zeros(nOB,nmBns);
    cov_Retr=zeros(nOB*(nOB-1)/2,nmBns);
    crr_Retr=zeros(nOB*(nOB-1)/2,nmBns); %padd any with 0 var with 0 corr
    denmCrr_Retr=zeros(nOB*(nOB-1)/2,nmBns);
    mn_spkCn=zeros(nOB,nmBns);
    vr_spkCn=zeros(nOB,nmBns);
    cov_spkCn=zeros(nOB*(nOB-1)/2,nmBns);
    crr_spkCn=zeros(nOB*(nOB-1)/2,nmBns); %padd any with 0 var with 0 corr
    denmCrr=zeros(nOB*(nOB-1)/2,nmBns);
    
    nTrial   = size(OBcount_allTrial_allCell,2); %all are 20
    
    
    for k=kDat
        % For each cell
        OBtemp = OBcount_allTrial_allCell(FirstSpon:LastEvok,:,k);
        OBtempHlf = OBcount_allTrial_allCellHlf(FirstSpon:LastEvok,:,k);
        
        mn_Orth(k,1:end-1)=(mean(OBtemp(1:end-1,1:10)')+ mean(OBtempHlf(1:end-1,1:10)')+ mean(OBtempHlf(2:end,1:10)'))/3;
        vr_Orth(k,1:end-1)=(var(OBtemp(1:end-1,1:10)')+var(OBtempHlf(1:end-1,1:10)')+var(OBtempHlf(2:end,1:10)'))/3;
        mn_Retr(k,1:end-1)=(mean(OBtemp(1:end-1,11:20)')+ mean(OBtempHlf(1:end-1,11:20)')+ mean(OBtempHlf(2:end,11:20)'))/3;
        vr_Retr(k,1:end-1)=(var(OBtemp(1:end-1,11:20)')+var(OBtempHlf(1:end-1,11:20)')+var(OBtempHlf(2:end,11:20)'))/3;
        mn_spkCn(k,1:end-1)=(mean(OBtemp(1:end-1,:)')+mean(OBtempHlf(1:end-1,:)')+mean(OBtempHlf(2:end,:)'))/3;
        vr_spkCn(k,1:end-1)=(var(OBtemp(1:end-1,:)')+var(OBtempHlf(1:end-1,:)')+var(OBtempHlf(1:end-1,:)'))/3;
        %last pt has only 2
        mn_Orth(k,end)=(mean(OBtemp(end,1:10))+ mean(OBtempHlf(end,1:10)))/2;
        vr_Orth(k,end)=(var(OBtemp(end,1:10))+var(OBtempHlf(end,1:10)))/2;
        mn_Retr(k,end)=(mean(OBtemp(end,11:20))+ mean(OBtempHlf(end,11:20)))/2;
        vr_Retr(k,end)=(var(OBtemp(end,11:20))+var(OBtempHlf(end,11:20)))/2;
        mn_spkCn(k,end)=(mean(OBtemp(end,:))+mean(OBtempHlf(end,:)))/2;
        vr_spkCn(k,end)=(var(OBtemp(end,:))+var(OBtempHlf(end,:)))/2;
        
        mean_sp1(k,:)=mean(OBtemp');
        mean_sp2(k,:)=mean(OBtempHlf');
        meanOrth1(k,:)=mean(OBtemp(:,1:10)');
        meanOrth2(k,:)=mean(OBtempHlf(:,1:10)');
        meanRetr1(k,:)=mean(OBtemp(:,11:20)');
        meanRetr2(k,:)=mean(OBtempHlf(:,11:20)');
    end
    
    %indicies for Cov (non-var)
    linInd=[];
    for rwI=1:(nOB-1)
        for clI=rwI+1:nOB
            linInd=[linInd; sub2ind([nOB nOB],rwI,clI)];
        end
    end
    
    %loop over time bins, not too inefficient
    for j=1:(nmBns-1)
        % 20 x nOB matrix
        tmpOB=squeeze(OBcount_allTrial_allCell(FirstSpon+(j-1),:,:)); %get all data
        tmpOB=tmpOB-repmat(mean_sp1(:,j)',nTrial,1); %centered
        covM_tmp=tmpOB'*tmpOB./(nTrial-1); %unbiased estim Cov; nOB x nOB matrix
        tmpVr_trun=vr_spkCn(:,j);
        tmpVr_trun(tmpVr_trun<1e-4)=1e-4; %lower Bound on Var for correl calc
            tmpOB2=squeeze(OBcount_allTrial_allCellHlf(FirstSpon+(j-1),:,:)); %get all data
            tmpOB2=tmpOB2-repmat(mean_sp2(:,j)',nTrial,1); %centered
            covM_tmp2=tmpOB2'*tmpOB2./(nTrial-1); %unbiased estim Cov; nOB x nOB matrix
        tmpOB3=squeeze(OBcount_allTrial_allCellHlf(FirstSpon+j,:,:)); %get all data
        tmpOB3=tmpOB3-repmat(mean_sp2(:,j+1)',nTrial,1); %centered
        covM_tmp3=tmpOB3'*tmpOB3./(nTrial-1); %unbiased estim Cov; nOB x nOB matrix
        corM_tmp=(covM_tmp+covM_tmp2+covM_tmp3)./(3*sqrt(tmpVr_trun*tmpVr_trun'));
        corM_tmpd=sqrt(tmpVr_trun*tmpVr_trun');
        
        cov_spkCn(:,j)=(covM_tmp(linInd)+covM_tmp2(linInd)+covM_tmp3(linInd))/3;
        crr_spkCn(:,j)=corM_tmp(linInd);
        denmCrr(:,j)=corM_tmpd(linInd);
        
        % Ortho: 10 x nOB matrix
        tmpOB=squeeze(OBcount_allTrial_allCell(FirstSpon+(j-1),1:10,:)); %get all data
        tmpOB=tmpOB-repmat(meanOrth1(:,j)',nTrial/2,1); %centered
        covM_tmp=tmpOB'*tmpOB./(nTrial/2-1); %unbiased estim Cov; nOB x nOB matrix
        tmpVr_trun=vr_Orth(:,j);
        tmpVr_trun(tmpVr_trun<1e-4)=1e-4; %lower Bound on Var for correl calc
            tmpOB2=squeeze(OBcount_allTrial_allCellHlf(FirstSpon+(j-1),1:10,:)); %get all data
            tmpOB2=tmpOB2-repmat(meanOrth2(:,j)',nTrial/2,1); %centered
            covM_tmp2=tmpOB2'*tmpOB2./(nTrial/2-1); %unbiased estim Cov; nOB x nOB matrix
        tmpOB3=squeeze(OBcount_allTrial_allCellHlf(FirstSpon+j,1:10,:)); %get all data
        tmpOB3=tmpOB3-repmat(meanOrth2(:,j+1)',nTrial/2,1); %centered
        covM_tmp3=tmpOB3'*tmpOB3./(nTrial/2-1); %unbiased estim Cov; nOB x nOB matrix
        corM_tmp=(covM_tmp+covM_tmp2+covM_tmp3)./(3*sqrt(tmpVr_trun*tmpVr_trun'));
        cov_Orth(:,j)=(covM_tmp(linInd)+covM_tmp2(linInd)+covM_tmp3(linInd))/3;
        crr_Orth(:,j)=corM_tmp(linInd);
        corM_tmpd=sqrt(tmpVr_trun*tmpVr_trun');
        denmCrr_Orth(:,j)=corM_tmpd(linInd);
        
        % Retro: 10 x nOB matrix
        tmpOB=squeeze(OBcount_allTrial_allCell(FirstSpon+(j-1),11:20,:)); %get all data
        tmpOB=tmpOB-repmat(meanRetr1(:,j)',nTrial/2,1); %centered
        covM_tmp=tmpOB'*tmpOB./(nTrial/2-1); %unbiased estim Cov; nOB x nOB matrix
        tmpVr_trun=vr_Retr(:,j);
        tmpVr_trun(tmpVr_trun<1e-4)=1e-4; %lower Bound on Var for correl calc
            tmpOB2=squeeze(OBcount_allTrial_allCellHlf(FirstSpon+(j-1),11:20,:)); %get all data
            tmpOB2=tmpOB2-repmat(meanRetr2(:,j)',nTrial/2,1); %centered
            covM_tmp2=tmpOB2'*tmpOB2./(nTrial/2-1); %unbiased estim Cov; nOB x nOB matrix
        tmpOB3=squeeze(OBcount_allTrial_allCellHlf(FirstSpon+j,11:20,:)); %get all data
        tmpOB3=tmpOB3-repmat(meanRetr2(:,j+1)',nTrial/2,1); %centered
        covM_tmp3=tmpOB3'*tmpOB3./(nTrial/2-1); %unbiased estim Cov; nOB x nOB matrix
        corM_tmp=(covM_tmp+covM_tmp2+covM_tmp3)./(3*sqrt(tmpVr_trun*tmpVr_trun'));
        cov_Retr(:,j)=(covM_tmp(linInd)+covM_tmp2(linInd)+covM_tmp3(linInd))/3;
        crr_Retr(:,j)=corM_tmp(linInd);
        corM_tmpd=sqrt(tmpVr_trun*tmpVr_trun');
        denmCrr_Retr(:,j)=corM_tmpd(linInd);
    end
    j=nmBns; %last bin; 2 pts
    tmpOB=squeeze(OBcount_allTrial_allCell(FirstSpon+(j-1),:,:)); %get all data
    tmpOB=tmpOB-repmat(mean_sp1(:,j)',nTrial,1); %centered
    covM_tmp=tmpOB'*tmpOB./(nTrial-1); %unbiased estim Cov; nOB x nOB matrix
    tmpVr_trun=vr_spkCn(:,j);
    tmpVr_trun(tmpVr_trun<1e-4)=1e-4; %lower Bound on Var for correl calc
    tmpOB2=squeeze(OBcount_allTrial_allCellHlf(FirstSpon+j-1,:,:)); %get all data
    tmpOB2=tmpOB2-repmat(mean_sp2(:,j)',nTrial,1); %centered
    covM_tmp2=tmpOB2'*tmpOB2./(nTrial-1); %unbiased estim Cov; nOB x nOB matrix
    corM_tmp=(covM_tmp+covM_tmp2)./(2*sqrt(tmpVr_trun*tmpVr_trun'));
    corM_tmpd=sqrt(tmpVr_trun*tmpVr_trun');
    cov_spkCn(:,j)=(covM_tmp(linInd)+covM_tmp2(linInd))/2;
    crr_spkCn(:,j)=corM_tmp(linInd);
    denmCrr(:,j)=corM_tmpd(linInd);
    
    % Ortho: 10 x nOB matrix
    tmpOB=squeeze(OBcount_allTrial_allCell(FirstSpon+(j-1),1:10,:)); %get all data
    tmpOB=tmpOB-repmat(meanOrth1(:,j)',nTrial/2,1); %centered
    covM_tmp=tmpOB'*tmpOB./(nTrial/2-1); %unbiased estim Cov; nOB x nOB matrix
    tmpVr_trun=vr_Orth(:,j);
    tmpVr_trun(tmpVr_trun<1e-4)=1e-4; %lower Bound on Var for correl calc
    tmpOB2=squeeze(OBcount_allTrial_allCellHlf(FirstSpon+j-1,1:10,:)); %get all data
    tmpOB2=tmpOB2-repmat(meanOrth2(:,j)',nTrial/2,1); %centered
    covM_tmp2=tmpOB2'*tmpOB2./(nTrial/2-1); %unbiased estim Cov; nOB x nOB matrix
    corM_tmp=(covM_tmp+covM_tmp2)./(2*sqrt(tmpVr_trun*tmpVr_trun'));
    cov_Orth(:,j)=(covM_tmp(linInd)+covM_tmp2(linInd))/2;
    crr_Orth(:,j)=corM_tmp(linInd);
    corM_tmpd=sqrt(tmpVr_trun*tmpVr_trun');
    denmCrr_Orth(:,j)=corM_tmpd(linInd);
    
    % Retro: 10 x nOB matrix
    tmpOB=squeeze(OBcount_allTrial_allCell(FirstSpon+(j-1),11:20,:)); %get all data
    tmpOB=tmpOB-repmat(meanRetr1(:,j)',nTrial/2,1); %centered
    covM_tmp=tmpOB'*tmpOB./(nTrial/2-1); %unbiased estim Cov; nOB x nOB matrix
    tmpVr_trun=vr_Retr(:,j);
    tmpVr_trun(tmpVr_trun<1e-4)=1e-4; %lower Bound on Var for correl calc
    tmpOB2=squeeze(OBcount_allTrial_allCellHlf(FirstSpon+j-1,11:20,:)); %get all data
    tmpOB2=tmpOB2-repmat(meanRetr2(:,j)',nTrial/2,1); %centered
    covM_tmp2=tmpOB2'*tmpOB2./(nTrial/2-1); %unbiased estim Cov; nOB x nOB matrix
    corM_tmp=(covM_tmp+covM_tmp2)./(2*sqrt(tmpVr_trun*tmpVr_trun'));
    cov_Retr(:,j)=(covM_tmp(linInd)+covM_tmp2(linInd))/2;
    crr_Retr(:,j)=corM_tmp(linInd);
    corM_tmpd=sqrt(tmpVr_trun*tmpVr_trun');
    denmCrr_Retr(:,j)=corM_tmpd(linInd);
end




