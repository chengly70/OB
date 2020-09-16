Scalar firing rate model using FI-curve from Li&Cleland model; same membrane time-constants, synaptic dynamics, etc.

Li, G. and T. A. Cleland (2013). A two-layer biophysical model of cholinergic neuromodulation in olfactory bulb. Journal of Neuroscience 33(7), 3037–3058. 
Li, G. and T. A. Cleland (2017). A coupled-oscillator model of olfactory bulb gamma oscillations. PLoS computational biology 13(11), e1005760.

call_mcRest.m -- example of a single calc for 1 point in parameter space (3D conductance and Time Vary Sigma), calls mc_Rest.m.
INPUT: indR index to specify what part of 3D conductance space, see lines 37-39 and dSeqHalton2.mat [which is random halton sequence sampling, see below] for how indR corresponds to a point in 3D conductance space

get_GsRand.m -- script detailing how random Gs are generated after getting good coup strengths from Monte Carlo of firing rate (mc_Rest.m) 
get_Ratios.m -- script detailing how random ratios are generated after getting ratios from Monte Carlo of firing rate

Main function to run firing rate model: mc_Rest.m
This function loads 2 mat-files for FI curves: dFI_coupldGc.mat and dFI_coupldPgc.mat (see /OB/ClelandMulticom/ for how created). 

INPUTS: Nc=# cells (7), Pstruct (struct of params; fields: tau_vec=membrane time-constants, CinMat=Nc x Nc background correlation matrix, cellTyp=Nc x 1 denoting 
cell type with 0=MC, 1=GC, 2=PGC), tmv=time-vector (-2:.001:2), muMat=Nc X length(tmv) for time-varying mean input, sig_vec1=Nc X 1 noise in Spont, sig_vec2=Nc X 1 noise in Evoked
N_relz=# of realizations, indR=index (1,2,..,10000) to denote region in 3D conduct space.

OUTPUTS: saves firing rate, variance and covariance (instantaneous and in 100ms windows) time-series in nameStr[indR].mat
See function comments for further details.

dSeqHalton2.mat is generated in Matlab via:

'>> hltob=haltonset(3,'Skip',1e3,'Leap',1e2);

'>> nvcSmp=hltob(1:10000,:);

'>> save dSeqHalton2 nvcSmp
