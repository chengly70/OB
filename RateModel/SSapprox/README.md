Contains steady-state approximation method, augmention of method in 
A. Barreiro & C. Ly, 2017. Practical approximation method for firing rate models of coupled neural networks with correlated inputs. 
Physical Review E (2), Vol. 96: pp 022413.
to allow larger coupling strengths

get_GsRand.m -- script detailing how random Gs are generated after getting good coup strengths from Monte Carlo of firing rate (mc_Rest.m) 

get_Ratios.m -- script detailing how random ratios are generated after getting ratios from Monte Carlo of firing rate

MAT Files: dGs_n.mat contains good coup strengths (Gs_good), each row has w_MG, w_GM, w_Gc. 
nRst242.mat is result of >> mc_Rest(242); has instance of results