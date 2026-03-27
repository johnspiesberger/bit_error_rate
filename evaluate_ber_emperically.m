% Evalutes bit error rate emperically.  Calls classify_symbols.m
%
% INPUTS
%
% fnam_sim_c_info_mat  String. Name of sim_c_info.mat file made and defined by simulate_c_information_measurement.m.
% 

% OUTPUTS
%
% ierr                 0: no error detected
%                      1; otherwise
%
function [ierr]=evaluate_ber_emperically(fnam_sim_c_info_mat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preliminary Steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
progname='evaluate_ber_emperically';
% fprintf(1,'%s: See jls7ub.mat to debug. pausing\n',progname);save jls7ub,pause




dotest=1;
if(dotest==1)
    fprintf(1,'%s: Test case. Press ENTER to continue\n',progname);pause
    
    fnam_sim_c_info_mat='sim_c_info.mat';

    t_classify=1.55e-05;
end

% Intialize outputs
ierr=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[a]=load(fnam_sim_c_info_mat);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Emperically derive BER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ber_emperical,ierr]=classify_symbols(a.h.r0,a.h.r1,a.h.dt,a.h.t_s,a.h.ind_direct_na_arrive,...
    a.h.ind_tau_alpha,a.h.ind_na_arrive_no_filters_no_interference,a.h.snr,a.invars.n_symbol_transmit, ...
    a.invars.n_samps_before_c,t_classify,a.h.t,[]);
if(ierr~=0),return,end


