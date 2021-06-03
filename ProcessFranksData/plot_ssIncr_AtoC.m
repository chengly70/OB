%script to plot pop-avg stats from ALL 8 mice

% stimt_clean - times of clean air stimuli (i.e. this is a control)
% stimt_eb - times of ethyl butyrate (our EB) stimuli 
% stimt_2hex - times of 2-hexanone stimuli 
% stimt_iso - times of isoamyl acetate stimuli 
% stimt_hexa - times of hexanal (our hexa) stimuli 
% stimt_et - times of ethyl tiglate stimuli 
% stimt_ea - times of ethyl acetate stimuli 

clrCd=[0 0 0];  %color of plots

flag_showEB=1; %show error bars (1) or not (0)
    stdf=0.05; %factor to mult std dev
flag_hlfSm=1; %half overlapping windows (1) or not (0)
nTrls=15; %same for all odors, all data
Twin=0.02; %20ms, get enough samples for p-values
sponTim=2.2; %have extra to cover 200ms before & after
evokTim=2; 
tme=(-sponTim+Twin:Twin:evokTim)';
LnmW=length(tme);

% specify which parts for spont & evoked (200ms before & after sniff)
iStrt=round((-.2-tme(1))/Twin)+1; %-0.2s
iDiv=round((0-tme(1))/Twin)+1; %0sec
iEv=round((.2-tme(1))/Twin)+1; %0.2s


for whch_odor=[1 4]  %0(no odor), 1=EB, 4=Hexanol
    
    % OUTPUTS, overwrite with each odor ---
    mnCntPop=[];
    vrCntPop=[];
    cvCntPop=[];
    nOBv=zeros(8,1);
    
    for whch_mouse=1:8
        switch whch_mouse
            case 1
                load('170608.mat')
            case 2
                load('170609.mat')
            case 3
                load('170613.mat')
            case 4
                load('170614.mat')
            case 5
                load('170618.mat')
            case 6
                load('170619.mat')
            case 7
                load('170621.mat')
            case 8
                load('170622.mat')
        end
        
        switch whch_odor
            case 0
                odStTimes=stimt_clean;
            case 1
                odStTimes=stimt_eb;
                clrCd=[0 0 0];  %color of plots
            case 2
                odStTimes=stimt_2hex;
            case 3
                odStTimes=stimt_iso;
            case 4
                odStTimes=stimt_hexa;
                clrCd=[0 0.5 0];  %color of plots
            case 5
                odStTimes=stimt_et;
            case 6
                odStTimes=stimt_ea;
        end
        
        
        
        numOB=length(obst);
        
        % TEMP outputs
        mnCnt=zeros(numOB,LnmW);
        vrCnt=zeros(numOB,LnmW);
        cvCnt=zeros(round(numOB*(numOB-1)/2),LnmW);
        
        Nct=cell(numOB,1);
        if(flag_hlfSm)
            Nht=cell(numOB,1); %count on half-grid pts
            up_nmTrl=3*nTrls;
            tmpM1=zeros(up_nmTrl,LnmW);
            tmpM2=zeros(up_nmTrl,LnmW);
        end
        
        for j=1:numOB
            tmp=obst{j};
            Nct{j}=zeros(nTrls,LnmW);
            if(flag_hlfSm)
                Nht{j}=zeros(nTrls,LnmW-1);
            end
            for k=1:nTrls
                id_get= (tmp>odStTimes(k)-sponTim)&(tmp<odStTimes(k)+evokTim);
                if(sum(id_get)>0)
                    spTims{k,j}=tmp(id_get)-odStTimes(k);
                else
                    spTims{k,j}=[];
                end
                
                Nct{j}(k,:)=histcounts(spTims{k,j},[tme-Twin;evokTim]);
                if(flag_hlfSm)
                    Nht{j}(k,:)=histcounts(spTims{k,j},tme-.5*Twin);
                end
            end
            if(flag_hlfSm)
                mnCnt(j,:)=mean( [Nct{j}; Nht{j} Nct{j}(:,end); Nct{j}(:,1) Nht{j}] );
                vrCnt(j,:)=var( [Nct{j}; Nht{j} Nct{j}(:,end); Nct{j}(:,1) Nht{j}] );
            else
                mnCnt(j,:)=mean(Nct{j});
                vrCnt(j,:)=var(Nct{j});
            end
        end
        
        %indicies for Cov (non-var)
        snm=1;
        %linInd=[];
        for rwI=1:(numOB-1)
            for clI=rwI+1:numOB
                %linInd=[linInd; sub2ind([numOB numOB],rwI,clI)];
                if(flag_hlfSm)
                    tmpM1=[Nct{rwI}; Nht{rwI} Nct{rwI}(:,end); Nct{rwI}(:,1) Nht{rwI}];
                    tmpM2=[Nct{clI}; Nht{clI} Nct{clI}(:,end); Nct{clI}(:,1) Nht{clI}];
                    
                    cvCnt(snm,:)=sum(tmpM1.*tmpM2)./(up_nmTrl-1)-up_nmTrl.*mnCnt(rwI,:).*mnCnt(clI,:)./(up_nmTrl-1);
                else
                    cvCnt(snm,:)=(sum(Nct{rwI}.*Nct{clI})-nTrls*mnCnt(rwI,:).*mnCnt(clI,:))./(nTrls-1);
                end
                
                snm=snm+1;
            end
        end
        
        %save data from each recording
        nOBv(whch_mouse)=numOB; %only when LOTS of mice, in order
        mnCntPop=[mnCntPop; mnCnt];
        vrCntPop=[vrCntPop; vrCnt];
        cvCntPop=[cvCntPop; cvCnt];
        
        
    end %for whch_mouse=1:8
    
    msp=mean(mnCntPop(:,iStrt:iDiv),2);
    mev=mean(mnCntPop(:,iDiv+1:iEv),2);
    vsp=mean(vrCntPop(:,iStrt:iDiv),2);
    vev=mean(vrCntPop(:,iDiv+1:iEv),2);
    csp=mean(cvCntPop(:,iStrt:iDiv),2);
    cev=mean(cvCntPop(:,iDiv+1:iEv),2);
    
    TotalOB = sum(nOBv); %get total # OB cells
    TotalPairs = size(cvCntPop,1); % get total # pairs
    
    Mean_PctIncr = sum(mev>msp)/TotalOB
    idCrs=mev>msp; %if ONLY restrict to cells that increase, but that's cheating
    %mean([msp mev]./Twin)
    
    Var_PctIncr = sum(vev>vsp)/TotalOB
    %mean([vsp vev])
    
    Cov_PctIncr = sum(cev>csp)/TotalPairs
    %mean([csp cev])
    
    % p-values reported in Table:
    [h,pValue_psth]=ttest2(msp,mev,'VarType','unequal'); %same with /Twin
    [h,pValue_Var]=ttest2(vsp,vev,'VarType','unequal');
    [h,pValue_Cov]=ttest2(csp,cev,'VarType','unequal');
    %for diff displays in command so user knows EB or HX
    if(whch_odor==1)  %EB
        pValue_psth_EB=pValue_psth
        pValue_Var_EB=pValue_Var
        pValue_Cov_EB=pValue_Cov
    elseif(whch_odor==4) %HX
        pValue_psth_HX=pValue_psth
        pValue_Var_HX=pValue_Var
        pValue_Cov_HX=pValue_Cov
    end
    
    %% box plots with mean+/-std
    figure
    hold on
    plot(mean([msp mev]./Twin),'s-','color',clrCd,'LineWidth',2,'MarkerSize',12,'MarkerFaceColor',clrCd)
    lw1=mean(msp./Twin)-stdf*std(msp./Twin); hi1=mean(msp./Twin)+stdf*std(msp./Twin);
    lw2=mean(mev./Twin)-stdf*std(mev./Twin); hi2=mean(mev./Twin)+stdf*std(mev./Twin);
    plot(ones(1,2),[lw1 hi1],'-','color',clrCd,'LineWidth',1)
    plot(2*ones(1,2),[lw2 hi2],'-','color',clrCd,'LineWidth',1)
    set(gca,'FontSize',24)
    box off
    ylabel('PSTH')
    axis([.8 2.2 8 12])
    
    figure
    hold on
    plot(mean([vsp vev]),'s-','color',clrCd,'LineWidth',2,'MarkerSize',12,'MarkerFaceColor',clrCd)
    lw1=mean(vsp)-stdf*std(vsp); hi1=mean(vsp)+stdf*std(vsp);
    lw2=mean(vev)-stdf*std(vev); hi2=mean(vev)+stdf*std(vev);
    plot(ones(1,2),[lw1 hi1],'-','color',clrCd,'LineWidth',1)
    plot(2*ones(1,2),[lw2 hi2],'-','color',clrCd,'LineWidth',1)
    set(gca,'FontSize',24)
    box off
    ylabel('Var')
    axis([.8 2.2 .14 .2])
    
    figure
    hold on
    plot(mean([csp cev]),'s-','color',clrCd,'LineWidth',2,'MarkerSize',12,'MarkerFaceColor',clrCd)
    lw1=mean(csp)-stdf*std(csp); hi1=mean(csp)+stdf*std(csp);
    lw2=mean(cev)-stdf*std(cev); hi2=mean(cev)+stdf*std(cev);
    plot(ones(1,2),[lw1 hi1],'-','color',clrCd,'LineWidth',1)
    plot(2*ones(1,2),[lw2 hi2],'-','color',clrCd,'LineWidth',1)
    set(gca,'FontSize',24)
    box off
    ylabel('Cov')
    axis([.8 2.2 5e-4 8e-3])
    
    if(whch_odor==4) %SAME size for both EB & HX
        TotalOB 
        TotalPairs
    end
    
end