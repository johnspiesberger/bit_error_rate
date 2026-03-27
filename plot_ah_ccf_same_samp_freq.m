% Takes in two signals sharing a sample frequency and plots the ccf and its abs hilbert transf, alongside the peak lag 
% time, t_pk_lag. t_pk_lag = lag time of peak of ccf's abs hilbert transf
% Note: this function is similar to tdoa_via_ccf_for_rec_pair_cst_c.m, but takes in pre-generated signals rather than
%       generating internally.
%
% INPUTS
%
% s1        1 x n.  s1(i) is pressure fluctuation of the first signal at a receiver.
% s2        1 x n.  s2(i)   "       "       "       "    second     "    the same receiver as s1
% fs        1 x 1.  the sampling frequency of the signals
% fign      1 x 1.  Figure number for plot. If empty, no plots shown
%
%
% OUTPUTS
%
% t_pk_lag  1 x 1.  the flipped peak lag time (s) between the two signals, where t_pk_lag > 0 inidicates s2 is delayed
%                   compared to s1, and t_pk_lag < 0 indicates s2 is advanced compared to s1.
%                   t_pk_lag = peak time of ccf's abs hilbert transf
% ierr      0: No error detected. 1: otherwise
%
function [t_pk_lag,ierr]=plot_ah_ccf_same_samp_freq(s1,s2,fs,fign)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preliminary Steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
progname='plot_ah_ccf_same_samp_freq.m';
%fprintf(1,'%s: See cags6122.mat to  debug. Pausing\n',progname);save cags6122;pause

dotest=0;
if(dotest==1)
    fprintf(1,'%s: Test case. Press Enter to continue\n',progname);pause
    
    c = 1500;
    fign_empty = [];
    scenario=1;
    switch scenario
        case 1
            src_xz=[-500 -70];
            rec_1_xz=[-50 -70];
            rec_2_xz=[50 -70];
            tau=0.025;
            omega=93*2*pi;
            g=5000;
            boundary_type=1;

            [h1,~]=exact_time_domain_solution_reflected_direct_paths(c,src_xz,rec_1_xz,tau,omega,g, ...
                boundary_type, fign_empty);
            [h2,~]=exact_time_domain_solution_reflected_direct_paths(c,src_xz,rec_2_xz,tau,omega,g, ...
                boundary_type, fign_empty);

            s1 = h1.s;
            s2 = h2.s;
            fs = g*omega/(2*pi);
           


        otherwise
        fprintf(1,'%s: software not written for scenario=%d\n',progname,scenario);
        ierr = 1;
        return
    end % switch scenario
    
    fign=[];
end %if(dotest==1)

% Initialize outputs
t_pk_lag = [];
ierr=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cross correlate the signals to get t_pk_lag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------------------------------------------------------------------------------------------------------
% ensure the signals are the same number of samples by zero-padding the end of the shorter
%-----------------------------------------------------------------------------------------------------------------------
%find signal lengths
len1 = length(s1);
len2 = length(s2);
len_diff = len1 - len2;
% zero pad the shorter signal. 
if len_diff > 0
    %len_diff > 0 means s2 is shorter.
    pad = zeros(1,len_diff);
    s2 = [s2 pad];
elseif len_diff < 0
    %len_diff < 0 means s1 is shorter.
    pad = zeros(1,-1*len_diff);
    s1 = [s1 pad];
end
%make sure signals are same length
if length(s1) ~= length(s2)
    fprintf(1, '%s: Warning. Signal lengths are different after padding. length(s1) = %d, while length(s2) = %d.\n', ...
        progname, length(s1),length(s2));
    ierr = 1;
    return
end
%-----------------------------------------------------------------------------------------------------------------------
% Compute cross-correlation between received signals and fix lag corresponding to max of abs value hilbert transform
%-----------------------------------------------------------------------------------------------------------------------
% Cross correlate s2 with s1 (other way around has sign flipped based on matlab xcorr)
[z]=xcorr(s2,s1,'normalized');
% Get abs value of hilbert transform of z
abs_h_z=abs(hilbert(z));
% Find maximum  peak
[~,ind_pk_max]=max(abs_h_z);
% get length of h_max.s)
ns=length(s1);
% Get lag axis for abs_h_z: lag_axis
lag_axis=1:length(abs_h_z);
% lag is zero at index ns of abs_h_z
lag_axis=(lag_axis - lag_axis(ns))*(1/fs);

% Derive lag of peak, switching the sign (t_m > 0 means signal arrives at receiver 2 later than receiver 1,
% corresponding to a lag z < 0)
t_pk_lag = lag_axis(ind_pk_max);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the ccf and abs hilbert transform if desired
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(~isempty(fign))
    figure(fign);
    clf;
    plot(lag_axis,z,'k-');
    hold on
    plot(lag_axis,abs_h_z,'g-');
    % Plot peak
    plot(lag_axis(ind_pk_max),abs_h_z(ind_pk_max),'r.');
    xlabel("lag (s)");
    ylabel('Abs value Hilbert Tranform CCF \color{red}  peak')
    title(['CCf of received waveforms. peak lag = ' num2str(t_pk_lag) ' seconds'])
    saveas(fign,'abs_hil_ccf','fig');
end %if(~isempty(fign))