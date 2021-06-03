%script to plot pop-avg TIME-VARYING stats from ALL 8 mice

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
Twin=0.1; %100ms for time-varying model match
sponTim=2.2;  %match anesthetized prep (much longer stim)
evokTim=2; 
tme=(-sponTim+Twin:Twin:evokTim)';
LnmW=length(tme);


for whch_odor=[1 4]  %0(no odor), 1=EB, 4=Hexanol; to compare with anesth
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
        ffCnt=zeros(numOB,LnmW);
        corrCnt=zeros(round(numOB*(numOB-1)/2),LnmW);
        denmS=zeros(round(numOB*(numOB-1)/2),LnmW);
        
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
    
    if(exist('h1'))
        figure(h1)
    else
        h1=figure;
    end
    hold on
    plot(tme,mean(mnCntPop./Twin),'color',clrCd,'LineWidth',2)
    if(flag_showEB)
        plot(tme,mean(mnCntPop./Twin)+stdf*std(mnCntPop./Twin),'--','color',clrCd,'LineWidth',1)
        plot(tme,mean(mnCntPop./Twin)-stdf*std(mnCntPop./Twin),'--','color',clrCd,'LineWidth',1)
        xtop=mean(mnCntPop./Twin)+stdf*std(mnCntPop./Twin);
        xbot=mean(mnCntPop./Twin)-stdf*std(mnCntPop./Twin);
        plot([tme(1) tme(1)],[xbot(1) xtop(1)],'--','color',clrCd,'LineWidth',1) %connect error bars
        plot([tme(end) tme(end)],[xbot(end) xtop(end)],'--','color',clrCd,'LineWidth',1) %connect error bars
    end
    set(gca,'FontSize',18)
    set(gca,'XLim',[-2 2])
    xlabel('Time (s)')
    ylabel('Firing Rate (Hz)')
    
    if(exist('h2'))
        figure(h2)
    else
        h2=figure;
    end
    hold on
    plot(tme,mean(vrCntPop),'color',clrCd,'LineWidth',2)
    if(flag_showEB)
        plot(tme,mean(vrCntPop)+stdf*std(vrCntPop),'--','color',clrCd,'LineWidth',1)
        plot(tme,mean(vrCntPop)-stdf*std(vrCntPop),'--','color',clrCd,'LineWidth',1)
        xtop=mean(vrCntPop)+stdf*std(vrCntPop);
        xbot=mean(vrCntPop)-stdf*std(vrCntPop);
        plot([tme(1) tme(1)],[xbot(1) xtop(1)],'--','color',clrCd,'LineWidth',1) %connect error bars
        plot([tme(end) tme(end)],[xbot(end) xtop(end)],'--','color',clrCd,'LineWidth',1) %connect error bars
    end
    set(gca,'FontSize',18)
    set(gca,'XLim',[-2 2])
    xlabel('Time (s)')
    ylabel('Var Spike Counts')
    %
    if(exist('h3'))
        figure(h3)
    else
        h3=figure;
    end
    hold on
    plot(tme,mean(cvCntPop),'color',clrCd,'LineWidth',2)
    if(flag_showEB)
        plot(tme,mean(cvCntPop)+stdf*std(cvCntPop),'--','color',clrCd,'LineWidth',1)
        plot(tme,mean(cvCntPop)-stdf*std(cvCntPop),'--','color',clrCd,'LineWidth',1)
        xtop=mean(cvCntPop)+stdf*std(cvCntPop);
        xbot=mean(cvCntPop)-stdf*std(cvCntPop);
        plot([tme(1) tme(1)],[xbot(1) xtop(1)],'--','color',clrCd,'LineWidth',1) %connect error bars
        plot([tme(end) tme(end)],[xbot(end) xtop(end)],'--','color',clrCd,'LineWidth',1) %connect error bars
    end
    set(gca,'FontSize',18)
    set(gca,'XLim',[-2 2])
    xlabel('Time (s)')
    ylabel('Cov Spike Counts')
    
end %for whch_odor=[1 4]