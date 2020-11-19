function [curve_corrected, coeff_corr] = deadtime_correction_matrix (Curves,NumChan,FirstCh,LastCh,AcqTime)
% deadtime_correction for optode g-SiPM 
% curve: matrix with dimensions n_ch x n_rep where: n_ch = number of channel; n_rep= repetition of the meas 
% NB: curve is the pure histogram, do not insert channel with flag/counter
% otherwise they will contribute to deadtime!!!!!!
% ch_first/ch_last : first/last channel of the TDC used
% T_meas_ms: measurement time (i.e. integration time) in ms
% LdS, 13/09/2019 --> SOLUS project --> optode
ns2sec = 1e-9;
NRep = size(Curves,3);
NDelays = size(Curves,1);
ch_usefull = LastCh - FirstCh + 1;
dead_time_coeff = ((NumChan - FirstCh):-1:(NumChan - LastCh)).* 0.65 + 30;
dead_time_matrix = repmat(dead_time_coeff,[NDelays 1 NRep]);
CW = squeeze(sum(Curves,2));
CW(CW<1)=1;
T_dead = squeeze(ceil(sum((dead_time_matrix.*Curves),2)))./squeeze(CW);
coeff_corr = repmat((1./(1-CW.*(T_dead*ns2sec/AcqTime))),[1 1 ch_usefull]);
coeff_corr = permute(coeff_corr,[1 3 2]);
curve_corrected = Curves .* coeff_corr;
