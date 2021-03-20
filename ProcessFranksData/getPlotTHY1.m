%script to process BoldenFranks data (Sci 2018), fix a specific intensity &
%           combine recorings
% for OB specifically
% Light pulses on for 1 sec, 10 trials for each light intensity.  
% Intensities increase from 3-6, 2-7, 2-9 (3 times), see Table2 pg 3 of crcns_pcx-1_data_descprtion.pdf
% 40+60+80+80+80 trials

cc=winter(10);

trls_use=zeros(5,10); %which trials to use to correspn to whch_Inten
nTrls=10;

Twin=0.1; %100ms
sponTim=2.2; 
evokTim=2; %stim stays on for 1 sec
tme=(-sponTim+Twin:Twin:evokTim)';
LnmW=length(tme);


for whch_Inten=4:6 %vary over 3 light-intensities used in Bolding&Franks 2018

trls_use=[ 10*(whch_Inten-3)+(1:10) ; repmat(10*(whch_Inten-2)+(1:10),4,1)];

flag_showEB=1; %show error bars (1) or not (0)
    stdf=0.05; %factor to mult std dev
flag_hlfSm=1; %half overlapping windows (1) or not (0)
flag_SS=0; %calc/show SS, FF & Correl calculations
flagSSbar=0; %shows time-avg lines in 2 stats; also calc t-test for stat sig

fExpData=cell(10,4); %know # rows & colmns in fid ahead of time

fid=fopen('ExperimentCatalog_THY1.txt');
tline = fgetl(fid); %skip first line
for k=1:10 
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    fExpData(k,:)=strsplit(tline);
    fExpData{k,1}=strrep(fExpData{k,1},'\','/'); %replace \ with / (Mac)
end
fclose(fid);

liOB=[]; %linear index in fExpData that has OB recordings
for j=1:10
    if(fExpData{j,4}=='B')
        liOB=[liOB; j];
    end
end

nmRecrd=length(liOB); %total number of recordings
nOBv=zeros(nmRecrd,1); %save each # of cells

% OUTPUTS, overwrite with each light intensity ---
mnCntPop=[];
vrCntPop=[];
cvCntPop=[];
FFpop=[]; %pop avg of each cell's FF
CRpop=[]; %pop avg of each pair's correlation
denmStd=[]; %sig_j*sig_k, for linear regr => Correl

for whichRecrd=liOB'

    if((whch_Inten==2 || whch_Inten==7) && whichRecrd==1) %no 2,7 intens in 1
        continue;
    end
    if((whch_Inten==8 || whch_Inten==9) && (whichRecrd==1 || whichRecrd==2)) %no 8 or 9 intens in 1,2
        continue;
    end
    
cFlName=[pwd,'/processed/',fExpData{whichRecrd,1}];

load([cFlName,'.efd'],'-mat'); %load efd structure: 
Stim_ts=efd.LaserTimes.LaserOn{1}; %get all times light turns on
clear efd;

load([cFlName,'.st'],'-mat'); %load SpikeTimes structure: tsec, units, Wave, stwarped

numOB=length(SpikeTimes.tsec)-1; 
spTims=cell(nTrls,numOB); %shift time so 0<->stim

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
    tmp=SpikeTimes.tsec{j+1}; %skip 1st cell, get all spk times from jth unit
    Nct{j}=zeros(nTrls,LnmW);
    if(flag_hlfSm)
        Nht{j}=zeros(nTrls,LnmW-1);
    end
    for k=1:nTrls 
        id_get= (tmp>Stim_ts(trls_use(whichRecrd,k))-sponTim)&(tmp<Stim_ts(trls_use(whichRecrd,k))+evokTim);
        if(sum(id_get)>0)
            spTims{k,j}=tmp(id_get)-Stim_ts(trls_use(whichRecrd,k));    
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

ffCnt=vrCnt./mnCnt; %fano factor
idbad=mean(ffCnt,2);
ffCnt=ffCnt(~isnan(idbad),:); %throwout entire cell/pair when denom=0

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
        corrCnt(snm,:)=cvCnt(snm,:)./sqrt(vrCnt(rwI,:).*vrCnt(clI,:));
        denmS(snm,:)=sqrt(vrCnt(rwI,:).*vrCnt(clI,:));
        
        snm=snm+1;
    end
end

idbad=mean(corrCnt,2);
corrCnt=corrCnt(~isnan(idbad),:); %throwout entire cell/pair when denom=0

clear SpikeTimes 

%save data from each recording
nOBv(whichRecrd)=numOB;
mnCntPop=[mnCntPop; mnCnt];
vrCntPop=[vrCntPop; vrCnt];
cvCntPop=[cvCntPop; cvCnt];
FFpop=[FFpop; ffCnt];
CRpop=[CRpop; corrCnt];
denmStd=[denmStd; denmS];

end %for whichRecrd=1:nmRecrd

if(flag_hlfSm)
    tmPlot=tme+Twin;
end

if(exist('h1'))
    figure(h1)
else
    h1=figure;
end

hold on
plot(tmPlot,mean(mnCntPop./Twin),'color',cc(whch_Inten,:),'LineWidth',2)
if(flag_showEB)
    plot(tmPlot,mean(mnCntPop./Twin)+stdf*std(mnCntPop./Twin),'--','color',cc(whch_Inten,:),'LineWidth',1)
    plot(tmPlot,mean(mnCntPop./Twin)-stdf*std(mnCntPop./Twin),'--','color',cc(whch_Inten,:),'LineWidth',1)
    xtop=mean(mnCntPop./Twin)+stdf*std(mnCntPop./Twin);
    xbot=mean(mnCntPop./Twin)-stdf*std(mnCntPop./Twin);
    plot([tmPlot(1) tmPlot(1)],[xbot(1) xtop(1)],'--','color',cc(whch_Inten,:),'LineWidth',1) %connect error bars
    plot([tmPlot(end) tmPlot(end)],[xbot(end) xtop(end)],'--','color',cc(whch_Inten,:),'LineWidth',1) %connect error bars
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
plot(tmPlot,mean(vrCntPop),'color',cc(whch_Inten,:),'LineWidth',2)
if(flag_showEB)
    plot(tmPlot,mean(vrCntPop)+stdf*std(vrCntPop),'--','color',cc(whch_Inten,:),'LineWidth',1)
    plot(tmPlot,mean(vrCntPop)-stdf*std(vrCntPop),'--','color',cc(whch_Inten,:),'LineWidth',1)
    xtop=mean(vrCntPop)+stdf*std(vrCntPop);
    xbot=mean(vrCntPop)-stdf*std(vrCntPop);
    plot([tmPlot(1) tmPlot(1)],[xbot(1) xtop(1)],'--','color',cc(whch_Inten,:),'LineWidth',1) %connect error bars
    plot([tmPlot(end) tmPlot(end)],[xbot(end) xtop(end)],'--','color',cc(whch_Inten,:),'LineWidth',1) %connect error bars
end
set(gca,'FontSize',18)
set(gca,'XLim',[-2 2])
xlabel('Time (s)')
ylabel('Var Spike Counts')

if(exist('h3'))
    figure(h3)
else
    h3=figure;
end
hold on
plot(tmPlot,mean(cvCntPop),'color',cc(whch_Inten,:),'LineWidth',2)
if(flag_showEB)
    plot(tmPlot,mean(cvCntPop)+stdf*std(cvCntPop),'--','color',cc(whch_Inten,:),'LineWidth',1)
    plot(tmPlot,mean(cvCntPop)-stdf*std(cvCntPop),'--','color',cc(whch_Inten,:),'LineWidth',1)
    xtop=mean(cvCntPop)+stdf*std(cvCntPop);
    xbot=mean(cvCntPop)-stdf*std(cvCntPop);
    plot([tmPlot(1) tmPlot(1)],[xbot(1) xtop(1)],'--','color',cc(whch_Inten,:),'LineWidth',1) %connect error bars
    plot([tmPlot(end) tmPlot(end)],[xbot(end) xtop(end)],'--','color',cc(whch_Inten,:),'LineWidth',1) %connect error bars
end
set(gca,'FontSize',18)
set(gca,'XLim',[-2 2])
xlabel('Time (s)')
ylabel('Cov Spike Counts')


end %end for whch_Inten=4:6, for varying light-intensity
