function []=mc_RateAwk(Nc,Pstruct,tmv,muMat,N_relz,indR,whichSet)
%[]=mc_RateAwk(Nc,Pstruct,tmv,muMat,N_relz,indR,whichSet)
%Simulate full network of Nc (input) WC, coupled/corrNoise, etc
%sig_vec, tau_vec are all Nc x 1 column vectors!
%F is calc from full model
%Gm is an NcxNc coupling matrix, with no autaptic (diag(Gm)=0)
%CinMat is an NcxNc correlation matrix; ones on diag and PSD!
% last 3 inputs: t, muMat is an Nc x Lt matrix, sig_vec=Nc x 1
% cells; t is Lt x 1, mu_t is an Nc X Lt matrix
% ASSUMING t=0 is when evoked, then noise becomes sig_vec2
%--- should check: (sig_vec,tau_vec,rv_vec,sp_vec1) are all Nc x 1
%---    CinMat is PSD, Gm is Nc x Nc ---
% OUTPUT: saves subset of results from cov_F,covF_w, mn_F
% similar to mc_RateHlf.m (no EB)  BUT using different indR (loading set here)

rng(indR*739+23) %seed random number generator
Twin=0.1; %in sec 

%synaptic parameters from LiCleland (1/sec)
alph_gaba=1000/1.25;
bet_gaba=1000/18;
alph_ampa=1000/1;
bet_ampa=1000/5.5;
alph_nmda=1000/52;
bet_nmda=1000/343;

% Other parameters
tau_vec     = Pstruct.tau_vec;
CinMat      = Pstruct.CinMat;
cellTyp     = Pstruct.cellTyp; %0=Excit (MC), 1=Inhib (GC), 2=Inhib (PGC)

%specify GM here, dependence on indR
%load dSeqHalton %gets nvcSmp [0,1] in 3 dims
load dSeqHalton2 %gets nvcSmp [0,1] in 3 dims, 10000
Gm=zeros(Nc,Nc);
gP2M=-5;
gM2P=5;

%run for subset (see calc_AwakeErrs.m for how to get dCovAwk.mat)
switch whichSet
    case 1 
        load('dCovAwk2.mat','indCT_sc') %has ALL indCT_n, indCT_sc, indCT_cs
        indL=indCT_sc(indR); %label of index, not 1:indR
        nameStr='sSm';%nameStr='sgSm';
        sig_vec1=0.25*[.75; 1; 1.5; 1.5; 1.5; 1; .75]; %spontaneous noise
        sig_vec2=sig_vec1;
    case 2
        load('dCovAwk2.mat','indCT_n')
        indL=indCT_n(indR); %label of index, not 1:indR
        nameStr='sVry'; %nameStr='sgVry'; 
        sig_vec1=0.25*[.75; 1; 1.5; 1.5; 1.5; 1; .75]; %spontaneous noise
        sig_vec2=sig_vec1;
        sig_vec2([1 2 6 7])=1.5*sig_vec2([1 2 6 7]); %evoked noise
    case 3
        load('dCovAwk2.mat','indCT_cs') 
        indL=indCT_cs(indR); %label of index, not 1:indR
        nameStr='sLrg'; %nameStr='sgLrg'; 
        sig_vec1=0.25*[.75; 1; 1.5; 1.5; 1.5; 1; .75]; %spontaneous noise
        sig_vec2=sig_vec1;
        sig_vec2([1 2 6 7])=1.5*sig_vec2([1 2 6 7]); %evoked noise
        sig_vec1=sig_vec2;
end


    gG2M=-5*2*nvcSmp(indL,1);
    gGcin=-2.5*2*nvcSmp(indL,2);
    gM2G=5*2*nvcSmp(indL,3);

gG2G=0; %not in Cleland, MarellaErmentrout, etc.
Gm=[0   gM2P 0     0 0 0 0; 
   gP2M 0    gG2M gGcin 0 0 0;
   0  gM2G 0   gG2G 0 0 0;
   0  gM2G   gG2G  0 gG2G gM2G 0;
   0 0 0       gG2G  0      gM2G 0;
   0 0 0      gGcin gG2M 0  gP2M;
   0 0 0 0 0           gM2P 0];
Gm(4,:)=.5*Gm(4,:); %common GC getting 2 copies of E/I inputs

%%!! assuming sig_vec constant!! otherwise, must change R_cor (noise) at all steps!
%sig_vec = sig_vec; 
% if(length(sig_vec1)~=Nc || length(sig_vec2)~=Nc || length(tau_vec)~=Nc || length(cellTyp)~=Nc ... 
%     || size(Gm,1)~=Nc || size(Gm,2)~=Nc || size(CinMat,1)~=Nc || size(muMat,1)~= Nc || size(muMat,2)~=length(tmv))
%     disp('Check params!');
%     return;
% end

dt=tmv(2)-tmv(1);
Lt=length(tmv);
sq_dt=1/sqrt(dt);

% load precomputed FI (transfer) curves
load('dFI_coupldGc.mat','frateMC','Ibg_v','frateGC');
input_F=Ibg_v./100; %scale down, extra intrins currents
FI_mg=frateMC;
FI_gc=frateGC;
load('dFI_coupldPgc.mat','frateMC'); %assuming SAME Ibg_v
%input_pg=Ibg_v; %assuming same as dFI_coupldGc
FI_mp=frateMC;

barR=max([FI_mg;FI_mp;FI_gc]); %max firing rate
%padd transfer/FI curves with 0s & large FR
tmp_f1=(-10:1:input_F(1)-1)';
input_F=[tmp_f1; input_F];
FI_mg=[zeros(size(tmp_f1)); FI_mg];
FI_gc=[zeros(size(tmp_f1)); FI_gc];
FI_mp=[zeros(size(tmp_f1)); FI_mp];
tmp_f1=(input_F(end)+2 : 2 : input_F(end)+10)';
input_F=[input_F; tmp_f1];
FI_mg=[FI_mg;FI_mg(end)*ones(size(tmp_f1))];
FI_gc=[FI_gc;FI_gc(end)*ones(size(tmp_f1))];
FI_mp=[FI_mp;FI_mp(end)*ones(size(tmp_f1))];

%EASIER/FASTEF when Ibg_v are the same (assuming it is)
F_Mat=zeros(length(input_F),Nc);
for j=1:Nc
    switch cellTyp(j)
        case 0 %MC cell
            F_Mat(:,j)=.5*FI_mp+.5*FI_mg;
        case 1 %GC cell
            F_Mat(:,j)=FI_gc;
        case 2 %PGC cell
            F_Mat(:,j)=FI_mp*exp(-.3);
    end
end

Tol=1e-5; %so syn filter doesn't go to -infinity
nmSt_g=round(-log(Tol)/(alph_gaba+bet_gaba)/dt); %t-steps back for GABA_A
nmSt_a=round(-log(Tol)/(alph_ampa+bet_ampa)/dt); %t-steps back for AMPA
nmSt_n=round(-log(Tol)/(alph_nmda+bet_nmda)/dt); %t-steps back for NMDA
maxSts=max([nmSt_a;nmSt_g;nmSt_n]);
tmv_g=.5*dt*(1:2:2*nmSt_g-1)'; %time vect used in syn filter
tmv_a=.5*dt*(1:2:2*nmSt_a-1)';
tmv_n=.5*dt*(1:2:2*nmSt_n-1)';
%following matrices used for syn; 
%  sum( M1_*r(j-1:-1:j-nmSt).*exp(M2_*r(j-1:-1:j-nmSt)) )
M1_g=sparse( alph_gaba/barR*diag(exp(-bet_gaba*tmv_g))*(diag([dt;.5*dt*ones(nmSt_g-1,1)])+diag([.5*dt*ones(nmSt_g-1,1)],-1)) );
M2_g=sparse(zeros(nmSt_g));
M2_g(3:end,2:end-1)=tril(ones(nmSt_g-2));
M2_g(:,1)=[.5;1.5*ones(nmSt_g-1,1)];
M2_g=-alph_gaba/barR*dt*M2_g;
    M1_a=sparse( alph_ampa/barR*diag(exp(-bet_ampa*tmv_a))*(diag([dt;.5*dt*ones(nmSt_a-1,1)])+diag([.5*dt*ones(nmSt_a-1,1)],-1)) );
    M2_a=sparse(zeros(nmSt_a));
    M2_a(3:end,2:end-1)=tril(ones(nmSt_a-2));
    M2_a(:,1)=[.5;1.5*ones(nmSt_a-1,1)];
    M2_a=-alph_ampa/barR*dt*M2_a;
M1_n=sparse( alph_nmda/barR*diag(exp(-bet_nmda*tmv_n))*(diag([dt;.5*dt*ones(nmSt_n-1,1)])+diag([.5*dt*ones(nmSt_n-1,1)],-1)) );
M2_n=sparse(zeros(nmSt_n));
M2_n(3:end,2:end-1)=tril(ones(nmSt_n-2));
M2_n(:,1)=[.5;1.5*ones(nmSt_n-1,1)];
M2_n=-alph_nmda/barR*dt*M2_n;

%temp to store prior firing rates
fr_back=zeros(N_relz,Nc,maxSts);

tauInv=1./tau_vec;
itauMat=repmat(tauInv',N_relz,1);
covM_xi1=CinMat.*(sig_vec1*sig_vec1'); %pure cov matrix
R_cor1=chol(covM_xi1); %R_cor'*randn(Nc,1) gives correct noise term (RHS)
% -- now for evoked sig_vec ---
covM_xi2=CinMat.*(sig_vec2*sig_vec2'); %pure cov matrix
R_cor2=chol(covM_xi2); %R_cor'*randn(Nc,1) gives correct noise term (RHS)
tId_evok=round((-tmv(1))/dt)+1; %index for evoked state, assuming at t=0

%Start sims
muT   = muMat(:,1); %use 1st time
xS=randn(N_relz,Nc)*diag(sig_vec1./sqrt(2*tau_vec))+repmat(muT',N_relz,1); %# of cells
FxM=zeros(N_relz,Nc);
synM=zeros(N_relz,Nc);
%do a t=.5 run to get correct IC starting point
for k=1:round(.5/dt)
    xi=randn(N_relz,Nc);
    
    %update synapses
    for j=1:Nc
        if(cellTyp(j)==0)%excit
            tmp_f1=squeeze(fr_back(:,j,1:nmSt_a)); %N_relz x nmSt
            tmp_f2=squeeze(fr_back(:,j,1:nmSt_n));
            synM(:,j)=sum( (M1_a*(tmp_f1')).*exp(M2_a*(tmp_f1')) )';
            synM(:,j)=synM(:,j)+...
                sum( (M1_n*(tmp_f2')).*exp(M2_n*(tmp_f2')) )';
        else %inhib; PGC or GC
            tmp_f1=squeeze(fr_back(:,j,1:nmSt_g));
            synM(:,j)=sum( (M1_g*(tmp_f1')).*exp(M2_g*(tmp_f1')) )';
        end
    end
    
    %eqns for cells 
    xS=xS+dt*itauMat.*( -xS+repmat(muT',N_relz,1)+sq_dt*(xi*R_cor1)+synM*(Gm') );
    
    %set current firing rate
    for j=1:Nc %FxM is N_relz x Nc
        FxM(:,j)=interp1(input_F,F_Mat(:,j),xS(:,j),'pchip');
    end
    %update prior firing rates
    fr_back(:,:,2:end)=fr_back(:,:,1:end-1);
    fr_back(:,:,1)=FxM;
end

Ltwin=floor((tmv(end)-tmv(1))/Twin); %# Twin windows
%twinv=(tmv(1)+Twin : Twin : tmv(end))';
nmSteps_tw=round(Twin/dt); %# dt steps in each Twin
hlfSteps=round(.5*nmSteps_tw); %# dt steps for half Twin (overlapping Twin)
% --- prior OUTPUTS. ---
mn_X=zeros(Nc,Lt);    %mean of Xj
cov_X=zeros(Nc,Nc,Lt);  %entire cov matrix of Xj (includes var)
mn_F=zeros(Nc,Lt);    %mean of F(Xj)
cov_F=zeros(Nc,Nc,Lt);  %entire cov matrix of Fj (includes var)
% -- repeat but for Twin --
covF_w=zeros(Nc,Nc,Ltwin);  %entire cov matrix of Fj (includes var)

%temp vars for Twin
tempF_prev=zeros(N_relz,Nc);
tempF_curr=zeros(N_relz,Nc);
tempFhlf_prev=zeros(N_relz,Nc);
tempFhlf_curr=zeros(N_relz,Nc);
ind_tw=1;

%in main loop, since N_relz x Nc, use (L*xi')'=xi*L', where L'=R_cor
% similarly for (Gm*FxMat')'=FxMat*Gm'

for k=1:Lt
    xi=randn(N_relz,Nc);
    if(k<tId_evok)
        R_cor=R_cor1;
    else
        R_cor=R_cor2;
    end
    %update synapses
    for j=1:Nc
        if(cellTyp(j)==0)%excit
            tmp_f1=squeeze(fr_back(:,j,1:nmSt_a)); %N_relz x nmSt
            tmp_f2=squeeze(fr_back(:,j,1:nmSt_n));
            synM(:,j)=sum( (M1_a*(tmp_f1')).*exp(M2_a*(tmp_f1')) )';
            synM(:,j)=synM(:,j)+...
                sum( (M1_n*(tmp_f2')).*exp(M2_n*(tmp_f2')) )';
        else %inhib; PGC or GC
            tmp_f1=squeeze(fr_back(:,j,1:nmSt_g));
            synM(:,j)=sum( (M1_g*(tmp_f1')).*exp(M2_g*(tmp_f1')) )';
        end
    end
    muT   = muMat(:,k);  %Use current time
    %eqns for cells 
    xS=xS+dt*itauMat.*( -xS+repmat(muT',N_relz,1)+sq_dt*(xi*R_cor)+synM*(Gm') );
    
    %set current firing rate
    for j=1:Nc %FxM is N_relz x Nc
        FxM(:,j)=interp1(input_F,F_Mat(:,j),xS(:,j),'pchip');
    end
    %update prior firing rates
    fr_back(:,:,2:end)=fr_back(:,:,1:end-1);
    fr_back(:,:,1)=FxM;
    
    %keep running sum of stats/dens
    mn_X(:,k)=mean(xS)';  
    mn_F(:,k)=mean(FxM)'; 
    cov_X(:,:,k)=xS'*xS./N_relz - mn_X(:,k)*mn_X(:,k)';
    cov_F(:,:,k)=FxM'*FxM./N_relz - mn_F(:,k)*mn_F(:,k)';
    
    % for Twin ..
    tempF_curr=tempF_curr+FxM;
    tempFhlf_curr=tempFhlf_curr+FxM;
    if( mod(k,nmSteps_tw)==0 ) %update prev; reset curr
        %shift & reset
        tempF_prev=tempF_curr;
        tempF_curr=zeros(N_relz,Nc);
        ind_tw=ind_tw+1;
    end
    if( mod(k,nmSteps_tw)==hlfSteps && ind_tw==1 ) %reset hlf b/c 1st .5Twin
        tempFhlf_curr=zeros(N_relz,Nc);
    elseif( mod(k,nmSteps_tw)==hlfSteps && ind_tw>1 ) %1.5Twin or larger
        % set cov's at previous ind_tw
        if(ind_tw==2)  %first one; only use hlfPrev & prev
            covF_w(:,:,1)=( tempF_prev'*tempF_prev./N_relz - mean(tempF_prev)'*mean(tempF_prev) + ...
                tempFhlf_prev'*tempFhlf_prev./N_relz - mean(tempFhlf_prev)'*mean(tempFhlf_prev) )./(2*nmSteps_tw^2);
        else %use all 3: hlfPrev,hlfCurr & prev
            
            covF_w(:,:,ind_tw-1)=( tempF_prev'*tempF_prev./N_relz - mean(tempF_prev)'*mean(tempF_prev) + ...
                tempFhlf_prev'*tempFhlf_prev./N_relz - mean(tempFhlf_prev)'*mean(tempFhlf_prev) + ...
                tempFhlf_curr'*tempFhlf_curr./N_relz - mean(tempFhlf_curr)'*mean(tempFhlf_curr) )./(3*nmSteps_tw^2);
        end
        tempFhlf_prev=tempFhlf_curr;
        tempFhlf_curr=zeros(N_relz,Nc);
    end
        
end
%update last Twin here; use hlfPrev & prev

covF_w(:,:,end)=( tempF_prev'*tempF_prev./N_relz - mean(tempF_prev)'*mean(tempF_prev) + ...
                tempFhlf_prev'*tempFhlf_prev./N_relz - mean(tempFhlf_prev)'*mean(tempFhlf_prev) )./(2*nmSteps_tw^2);

%only save subset of results, for memory purposes
mcFm=.5*mn_F(2,:)+.5*mn_F(6,:);
varFRt=.5*squeeze(cov_F(2,2,:))' + 5*squeeze(cov_F(6,6,:))';
varFRt_tw=.5*squeeze(covF_w(2,2,:))' + .5*squeeze(covF_w(6,6,:))';
covFRt=squeeze(cov_F(2,6,:))'; %exact same as (6,2), x'*x
covFRt_tw=squeeze(covF_w(2,6,:))';
covFRt_tw(1)=2*covFRt_tw(1); %scale to match other bins
varFRt_tw(1)=2*varFRt_tw(:,1);

save([pwd,nameStr,num2str(indR)],'mcFm','varFRt','varFRt_tw','covFRt','covFRt_tw')

