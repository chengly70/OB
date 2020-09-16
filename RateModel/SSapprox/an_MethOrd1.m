function [convged,Corr_valid,cov_Fa,mn_Fa,cov_Xa,mn_Xa,mean_all]=an_MethOrd1(Nc,Pstruct,mnXan,mu_vec,sig_vec)
%[convged,Corr_valid,cov_Fa,mn_Fa,cov_Xa,mn_Xa,mean_all]=an_MethOrd1(Nc,Pstruct,mnXan,mu_vec,sig_vec,CinMat)
%sig_vec, tau_vec are all Nc x 1 column vectors!
%F is calc from full model
%Gm is an NcxNc coupling matrix, with no autaptic (diag(Gm)=0)
%CinMat is an NcxNc correlation matrix; ones on diag and PSD!
% last 3 inputs: t, muMat is an Nc x Lt matrix, sig_vec=Nc x 1
% cells; t is Lt x 1, mu_t is an Nc X Lt matrix
% ASSUMING t=0 is when evoked, then noise becomes sig_vec2
%OUTPUTS: convged: indicates whether self-consis system converged (1) or not (0)
% Corr_valid: indicates is corr-mat of actvitiy (X) is positive def (1) or not (0)
% cov_Fa = Nc x Nc covariance matrix of firing rate activity
% mn_Fa = Nc x 1 vector of mean firing rates
% cov_Xa = Nc x Nc covariance matrix of activity [Xj]
% mn_Xa = Nc x 1 vector of mean activity [Xj]
% mean_all = Nc x ? matrix to show mn_Xa at EACH step of iteration
% extention of an_Method.m, looking at deviations around mnXan[Sp/Ev]

Num_it=50; %(total) MAX # of iterations
Tolnc=1e-4; %tolerance level for convergence

%synaptic parameters from LiCleland (1/sec)
alph_gaba=1000/1.25;
bet_gaba=1000/18;
alph_ampa=1000/1;
bet_ampa=1000/5.5;
alph_nmda=1000/52;
bet_nmda=1000/343;

% Other parameters
tau_vec     = Pstruct.tau_vec;
Gm          = Pstruct.Gm;
CinMat      = Pstruct.CinMat;
cellTyp     = Pstruct.cellTyp; %0=Excit (MC), 1=Inhib (GC), 2=Inhib (PGC)

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

% --- 

id1=[];
id2=[];
for j=1:(Nc-1)
    id1=[id1; (1:j)'];
    id2=[id2; (j+1)*ones(j,1)];
end
ind_UpTri=sub2ind([Nc Nc],id1,id2); %indices upper triang Corr/Cov matrix

sclSig=1/sqrt(2);
%need discret norm pdfs to do analytic calcs
d_nrm=0.01;
norm_dom=(-3:d_nrm:3)';
len_nrm=length(norm_dom);
normDom_mat=repmat(norm_dom,1,Nc); %matrx len_nrm x Nc, use in F_mn/vr calcs
norm01_pdf=normpdf(norm_dom,0,1); %stand. normal
nrm_varF=repmat(norm_dom,len_nrm,1); %[x1(1),x1(2),...]
nrm_varS=reshape(repmat(norm_dom',len_nrm,1),len_nrm*len_nrm,1); %[x2(1),x2(1),...]    
%numerical correl 2D gaussian
crMat_p=triu(CinMat-diag(diag(CinMat))); %entries required for correl 2D Gaussians
ind_nnzCr=find(crMat_p);   %linear indices of non-zero correl
[rw_ind,cl_ind]=ind2sub([Nc Nc],ind_nnzCr); %row and column indices of non-zero correl
num_Crs=length(ind_nnzCr); %number of non-zero correl
ind_zrCr=find(CinMat==0); %linear indices of Zero correl
% [rz_i,cz_i]=ind2sub([Nc Nc],ind_zrCr); %row and column indices of Zero correl
% zrInd_keep=rz_i<cz_i;   %only keep upper triangular
% rz_i=rz_i(zrInd_keep);  %only keep upper triangular
% cz_i=cz_i(zrInd_keep);  %only keep upper triangular
% num_Zrcr=length(rz_i); 
if(num_Crs~=0) %only eval if have background correl
    norm2D_all=zeros(Nc,Nc,len_nrm*len_nrm); %3D matrix, all correl2Dgaussians
    for j=1:num_Crs
        sigTmp=[1 CinMat(rw_ind(j),cl_ind(j));  CinMat(rw_ind(j),cl_ind(j)) 1];
        norm2D_all(rw_ind(j),cl_ind(j),:)=mvnpdf([nrm_varF nrm_varS],[0 0],sigTmp); %0 mean, unit var, cov=CinMat(j,k)
    end
end
norm2D_sing=zeros(len_nrm*len_nrm,1);
tauInv_scl=1./(repmat(tau_vec,1,Nc)+repmat(tau_vec',Nc,1)); %scale cov_Xa by this, Nc x Nc

%-- temp variables --
Finst_1D=zeros(len_nrm,Nc);
Hinst_1D=zeros(len_nrm,Nc);
synM=zeros(Nc,1); 

%OUTPUTS; iterative solver, start at uncoupled values
mn_Xa=mnXan;    %mean of Xj, start at sigma=0
cov_Xa=tauInv_scl.*diag(sig_vec)*CinMat*diag(sig_vec);  %entire cov matrix of Xj (includes var)
corr_Xa=zeros(Nc,Nc); %entire correl matrix of Xj (ones on diagonal)
mn_Fa=zeros(Nc,1);    %mean of F(Xj)
cov_Fa=zeros(Nc,Nc);  %entire cov matrix of Fj (includes var)
%corr_Fa=zeros(Nc,Nc); %entire correl matrix of Fj (ones on diagonal)
%simulation does not converge by default
convged=(1==0);
Corr_valid=(1==1); %logical to indicate whether correl was inbounds (1) or not (0)

%save all of the means at each step
mean_all=mn_Xa;
mean_tmp=zeros(Nc,1);
%only save previous var/cov to check convergence
var_prev=diag(cov_Xa);
var_tmp=zeros(Nc,1);
cov_prev=nonzeros(triu(cov_Xa)-diag(var_tmp));
cov_tmpT=zeros(Nc*(Nc-1)*.5,1); %check covX (non-var) values
Mfsq=zeros(Nc,Nc); %matrix for cov of coupling direct terms
Mnf=zeros(Nc,Nc);  %matrix for cov of coupling with noise term

for j=2:Num_it
    %instance of F(sig(j-1)*y+mu(j-1)), len_nrm x Nc matrix for 1D-integ
    for k=1:Nc
        Finst_1D(:,k)=interp1(input_F,F_Mat(:,k),norm_dom*sqrt(var_prev(k))+repmat(mean_all(k,j-1)-mnXan(k),len_nrm,1),'pchip');
        %update synapses
        if(cellTyp(k)==0)%excit
            Hinst_1D(:,k)=alph_ampa/barR*Finst_1D(:,k)./(bet_ampa+alph_ampa/barR*Finst_1D(:,k)) + ...
                alph_nmda/barR*Finst_1D(:,k)./(bet_nmda+alph_nmda/barR*Finst_1D(:,k));
            synM(k)=norm01_pdf'*( Hinst_1D(:,k) )*d_nrm;
        else %inhib; PGC or GC
            Hinst_1D(:,k)= alph_gaba/barR*Finst_1D(:,k)./(bet_gaba+alph_gaba/barR*Finst_1D(:,k));
            synM(k)=norm01_pdf'*( Hinst_1D(:,k) )*d_nrm;
        end
    end
    
    %firing rate approx
%     F_mn=norm01_pdf'*Finst_1D*d_nrm; 
%     F_mn=F_mn'; %1 x Nc before this command
    F_secm=norm01_pdf'*( Hinst_1D.^2 )*d_nrm;
    F_secm=F_secm'; %1 x Nc before this command; second moment of F
        
    mean_tmp=mu_vec-mnXan+Gm*synM;
    
    %___ form the 2 large Nc x Nc matrices, using var_prev,mean_all(:,j-1) ___
    Mfsq=diag(F_secm-synM.^2); %set the diagonalS
    Mnf=diag(sclSig*norm01_pdf'*(Hinst_1D.*(normDom_mat*diag(sig_vec)))*d_nrm);
    %Mfsq and Mnf have save sparsity structure as CinMat
    %this for-loop will only eval if have background correl (CinMat not diag)
    for k=1:num_Crs 
        Mfsq(rw_ind(k),cl_ind(k))=...
            sum(repmat(Hinst_1D(:,rw_ind(k))-synM(rw_ind(k)),len_nrm,1).*reshape(repmat(Hinst_1D(:,cl_ind(k))'-synM(cl_ind(k)),len_nrm,1),len_nrm*len_nrm,1)...
            .*squeeze(norm2D_all(rw_ind(k),cl_ind(k),:)))*d_nrm*d_nrm;
        Mfsq(cl_ind(k),rw_ind(k))=Mfsq(rw_ind(k),cl_ind(k)); %make symmetric, this is faster then outside for-loop
        
        Mnf(rw_ind(k),cl_ind(k))=sclSig*sig_vec(cl_ind(k))*...
            sum(repmat(Hinst_1D(:,rw_ind(k)),len_nrm,1).*nrm_varS...
            .*squeeze(norm2D_all(rw_ind(k),cl_ind(k),:)))*d_nrm*d_nrm;
        %set the transpose entry too
        Mnf(cl_ind(k),rw_ind(k))=sclSig*sig_vec(rw_ind(k))*...
            sum(repmat(Hinst_1D(:,cl_ind(k)),len_nrm,1).*nrm_varS...
            .*squeeze(norm2D_all(rw_ind(k),cl_ind(k),:)))*d_nrm*d_nrm;
    end
    
    %perform the calculation
    cov_Xa=tauInv_scl.*(diag(sig_vec)*CinMat*diag(sig_vec)+Gm*Mnf+Mnf'*Gm'+Gm*Mfsq*Gm');
    
    %set tmp variables and updates
    mean_all=[mean_all mean_tmp]; %dynamic increase; can overwrite if want
    var_tmp=diag(cov_Xa);
    cov_tmpT=triu(cov_Xa)-diag(var_tmp);
    cov_tmpT=cov_tmpT(ind_UpTri);
    
    %check if iteration converges
    if(j>4 && norm(mean_all(:,j-1)-mean_all(:,j))<Tolnc && norm(var_prev-var_tmp)<Tolnc && norm(cov_tmpT-cov_prev)<2*Nc*Tolnc)
        convged=~convged;
        break; %all quantities have converged
    end 
    
    %Update prev variables, storing previous result
    var_prev=var_tmp; 
    cov_prev=cov_tmpT;
end

if(convged==0)
    
    return; %leave function; has not converged
    
else %calc rest of the stats
    mn_Xa=mean_all(:,end) + mnXan; %add on baseline mnXan
    varX_matr=var_tmp*var_tmp'; %used to divide to get correl matrix
    %calc the corr(Xj,Xk)
    corr_Xa=cov_Xa./sqrt(varX_matr);
        
    %check if correl matrix is positive definite (assuming semi-def is unlikely)
    [L,p_chol]=chol(corr_Xa);
    if(p_chol>0)
        Corr_valid=~Corr_valid;
        return; %exit function; corr structure is invalid
    end
    
    %calc firing stats, assuming (X_j,Y_k) are multivariate Normal
    for k=1:Nc
        Finst_1D(:,k)=interp1(input_F,F_Mat(:,k),norm_dom*sqrt(2*var_tmp(k))+repmat(mn_Xa(k),len_nrm,1),'pchip');
    end
    mn_Fa=norm01_pdf'*Finst_1D*d_nrm; %mean firing rate
    mn_Fa=mn_Fa'; %make it Nc x 1
    
    %set the variance of the firing rate
    cov_Fa=diag( (norm01_pdf'*( Finst_1D-repmat(mn_Fa',len_nrm,1) ).^2*d_nrm)' );
    
    %Recalc numeric 2D-Gauss with FINAL cov/corr, calc all .5Nc(Nc-1) of them
    %not saving 2D-Gauss; finish calculation here, cov_Fa off-diagonal
    for colInd=2:Nc
        for rowInd=1:(colInd-1)
            sigTmp=[1 corr_Xa(rowInd,colInd); corr_Xa(rowInd,colInd)  1];
            norm2D_sing=mvnpdf([nrm_varF nrm_varS],[0 0],sigTmp); %shift mean, unit var, cov=corr_Xa

            cov_Fa(rowInd,colInd)=...
            sum(repmat(Finst_1D(:,rowInd)-mn_Fa(rowInd),len_nrm,1).*reshape(repmat(Finst_1D(:,colInd)'-mn_Fa(colInd),len_nrm,1),len_nrm*len_nrm,1)...
            .*norm2D_sing)*d_nrm*d_nrm;
        
            cov_Fa(colInd,rowInd)=cov_Fa(rowInd,colInd); %make symmetric
        end
    end
end

