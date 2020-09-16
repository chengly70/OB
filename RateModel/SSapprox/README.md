Contains steady-state approximation method, augmention of method in 
A. Barreiro & C. Ly, 2017. Practical approximation method for firing rate models of coupled neural networks with correlated inputs. 
Physical Review E (2), Vol. 96: pp 022413.
to allow larger coupling strengths

call_AnMeth.m -- script to call SS approximation method, calls anMeanOnly.m and an_MethOrd1.m

anMeanOnly.m -- function to iteratively solve for shifted mean with no noise

an_MethOrd1.m -- derived from original method in Barreiro & Ly '17 but augmented to include specific OB firing rate model

get_GsRand.m -- script detailing how random Gs are generated after getting good coup strengths from Monte Carlo of firing rate (mc_Rest.m) 

get_Ratios.m -- script detailing how random ratios are generated after getting ratios from Monte Carlo of firing rate

MAT Files: dGs_n.mat contains good coup strengths (Gs_good), each row has w_MG, w_Gc, w_GM. 
nRst242.mat is result of >> mc_Rest(242); has one instance of results
