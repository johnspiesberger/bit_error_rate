% Run exact_time_domain_solution_reflected_direct_paths.m for one source and two receivers, then cross correlate the
% time series received by each to get the TDoA
% Note: this function is similar to plot_ah_ccf_same_samp_freq.m, but generates signals internally rather than requiring
%       pregenerated signals
%
% INPUTS
%
% c               1 x 1.  in-situ speed of sound in water (m/s)
% src_xz          1 x 1.  Source is located at  horizontal coordinate x=0 and  vertical cooordinate, z=src_z (m).
%                         z is 0 at ocean's surface and is positive up.  source also located at y=0,  in 3D space.
% rec_1_xz        1 x 2.  The first receiver's horizontal (rec_1_xz(1)) and vertical coordinate (rec_1_xz(2)). Same 
%                         coordinate system as for src_z             
% rec_2_xz        1 x 2.  Same arrangment as above, but for the second receiver
% tau             1 x 1.  emitted pulse is exp(-(t/tau)^2)*cos(omega*t) with tau in seconds and t in seconds
% omega           1 x 1.  See tau above. Units are rad/s
% samp_freq_mult  1 x 1.  Sample frequency of synthsized time series will equal samp_freq_mult * omega/(2*pi).
%                         Must be >=2.
% boundary_type   1 x 1.  1: boundary is like surface of ocean where there are no pressure fluctuations.
%                            There is a direct path and a path reflecting from surface 
%                         2: boundary is hard. Reflected path bounces from hard flat rocks on bottom of ocean.
% fign            1 x 1.  Figure number for plot. If empty, no plots shown
%
% OUTPUTS
%
% dt              1 x 1. The tdoa of the signal emitted by the source to the receivers (s). Equal to the peak
%                        lag produced by the cross correlation of the two received signals
% h1              Structured array containing data for receiver 1, with fields:
%
%   -taxis        1 x n.  taxis(i) is i'th time of data at receiver (s)
%   -s_emit       1 x n.  s_emit(i) is emitted signal at source at time  t(i).  This is the signal from the source
%                         in absence of any other signals.
%   -s            1 x n.  s(i) is pressure fluctuation of signal at receiver. s(i) is at  time taxis(i).
%   -c3d          1 x 1. 3D effective speed (m/s) measured from time of emitted  peak
%                        to time of  biggest received peak. Equals l_d/(t_r-t_emit)
%   -t_emit       1 x 1. Time (s) when emitted signal reaches peak
%   -t_r          1 x 1. Time (s) when received signal reaches maximum
%   -l_direct     1 x 1. Pythagorean distance (m)  between source and receiver.
% h2              Same as h1, except for receiver 2
% ierr            0: No error detected. 1: otherwise
%
function [dt,h1, h2,ierr]=tdoa_via_ccf_for_rec_pair_cst_c(c,src_z,rec_1_xz,rec_2_xz,tau,omega,samp_freq_mult, ...
    boundary_type,fign)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preliminary Steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
progname='tdoa_via_ccf_for_rec_pair_cst_c.m';
% fprintf(1,'%s: See cags6111.mat to  debug. Pausing\n',progname);save  cags6111;pause

dotest=1;
if(dotest==1)
    fprintf(1,'%s: Test case. Press Enter to continue\n',progname);pause
    c=1480;
    
    scenario=1;
    switch scenario

        case 1
            src_z=-20;
            rec_1_xz=[15.3 -15.7];
            rec_2_xz=[3 -8];
            tau=0.025;
            omega=100*2*pi;
            samp_freq_mult=10000;
            boundary_type=1;
        case 2
            src_z=-20;
            rec_1_xz=[15.3 -15.7];
            rec_2_xz=[3 -8];
            tau=0.025;
            omega=100*2*pi;
            samp_freq_mult=10000;
            boundary_type=1;

        otherwise
        fprintf(1,'%s: software not written for scenario=%d\n',progname,scenario);
        pause
    end % switch scenario
    
    fign=1;
end

% Initialize outputs
dt = [];
h1.taxis=[];
h1.s_emit=[];
h1.s=[];
h1.c3d=[];
h1.t_emit=[];
h1.t_r=[];
h1.l_direct=[];
h2 = h1;
ierr=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get dt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------------------------------------------------------------------------------------------------------
% simulate the signals
%-----------------------------------------------------------------------------------------------------------------------
% decide if plots from exact_time_domain_solution_reflected_direct_paths.m are wanted using fign1 and fign2
fign1 = [];
fign2 = [];
%simulate signals
g = samp_freq_mult;
[h1,ierr1]=exact_time_domain_solution_reflected_direct_paths(c,src_z,rec_1_xz,tau,omega,g,boundary_type,fign1);
if ierr1 == 1
    fprintf(1, '%s: exact_time_domain_solution_reflected_direct_paths.m returns ierr = 1 for\n',progname)
    fprintf(1, '    source at [0m, %fm] and receiver 1 at [%fm,%fm].\n',src_z, rec_1_xz(1), rec_1_xz(2))
    ierr = 1;
    return 
end     % ierr1 == 1

[h2,ierr2]=exact_time_domain_solution_reflected_direct_paths(c,src_z,rec_2_xz,tau,omega,g,boundary_type,fign2);
if ierr2 == 1
    fprintf(1, '%s: exact_time_domain_solution_reflected_direct_paths.m returns ierr = 1 for\n',progname)
    fprintf(1, '    source at [0m, %fm] and receiver 1 at [%fm,%fm].\n',src_z, rec_2_xz(1), rec_2_xz(2))
    ierr = 1;
    return 
end     %ierr2 == 1

%-----------------------------------------------------------------------------------------------------------------------
% Compute cross-correlation between received signals and fix lag corresponding to max of abs value hilbert transform
%-----------------------------------------------------------------------------------------------------------------------
fprintf(1, '%s: WARNING. Check whether lag should be positive or negative, as MATLAB xcorr may have it flipped.', ...
    progname)
% Cross correlate
[z]=xcorr(h1.s,h2.s,'none');
% [z]=xcorr(h_max.s_emit,h_max.s_emit,'biased');
% Get abs value of hilbert transform of z
abs_h_z=abs(hilbert(z));
% Find maximum  peak
[~,ind_pk_max]=max(abs_h_z);
% get length of h_max.s)
ns=length(h1.s);
% Get lag axis for abs_h_z: lag_axis
lag_axis=1:length(abs_h_z);
% lag is zero at index ns of abs_h_z
lag_axis=(lag_axis - lag_axis(ns))*diff(h1.taxis(1:2));
% Derive lag of peak
lag_pk_max=lag_axis(ind_pk_max);


% tdoa is the peak lag time
dt = lag_pk_max;
fprintf(1,'Using cross-correlation of waveforms received by receiver pair, TDoA = peak  of Hilbert transf= %f s\n', ...
    dt);

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
    plot(lag_pk_max,abs_h_z(ind_pk_max),'r.');
    xlabel("lag (s)");
    ylabel('Abs value Hilbert Tranform CCF \color{red}  peak')
    title(['CCf of received waveforms. TDoA = ' num2str(dt) ' seconds'])
    saveas(fign,'abs_hil_ccf','fig');
end %if(~isempty(fign))