% Simulates experiment where the speed of information is derived from measurements. Follows ideas from
% Stenner, Gauthier, and Neifeld letters to nature, vol 425, p. 695-8, 2003
%
% Program calls these subroutines
%
% design_bandpass_filter.m
% filt_the_data.m
% simulate_two_waveforms.m
%
% INPUTS
%
% invars            Structured array with fields:
%
% f                 Defined by simulate_two_waveforms.m
% tau               "
% tau_after         "
% amp_big           "
% f_samp            "
% n_tau_to_zero     "
%
% stopband_freq1    Defined by design_bandpass_filter.m
% passband_freq1    "
% passband_freq2    "
% stopband_freq2    "
% stopband_atten1   "
% passband_ripple   "
% stopband_atten2   "
% design_method     "
% fign_filter       1 x 1.  Specifies matlab figure number to show bandpass filter. If empty, plot not shown.
%
% receiver_filter   Structured array. Fields are same as above except pertain to the filter at the receiver.
%   -stopband_freq1    
%   -passband_freq1    
%   -passband_freq2    
%   -stopband_freq2    "
%   -stopband_atten1   "
%   -passband_ripple   "
%   -stopband_atten2   "
%   -design_method     "
%   -fign_filter       1 x 1.  Specifies matlab figure number to show bandpass filter. If empty, plot not shown.
%
% src_xz            1 x 2. (x,z) coordinate of source (m). z positive up and zero at boundary
% rec_xz            1 x 2. "                   receiver (m). "
% boundary_type     1 x 1. 1: zero pressure perturbation at surface
%                          2: No normal velocity at boundary (hard boundary)
% c                 1 x 1. Phase speed of signal (m/s).
% method_interp     String.  Equals method for interp1.m. E.g. 'linear'
% energy_loss       1 x 1.
%                   1: Energy decreases according to spherical spreading of energy with distance of path from source
%                   2: No energy decrease with distance.  This is what happens with a laser before the laser beam
%                      starts to spread out significantly
% snr               1 x 1. Signal-to-noise ratio (dB) at receiver for symbol 1 at its peak value.
% n_symbol_transmit 1 x 1.  Number of times symbol 0 and symbol 1 are transmitted
%                   (each transmitted n_symbol_transmit times). Need enough of these to
%                   get a reliable bit error rate (BER). e.g. 100
% n_samps_before_c
%                   1 x 1. Defined by bit_error_rate.m.
% seedd             1 x 1. Non-negative integer.  Used to initialize random number generators
% ber_thresh        1 x 1. Bit error rate for accepting significant classification of symbols. paper 1 uses 0.1
%
% OUTPUTS
%
% ber_info_interfere
%                 Structured array made and defined by ber_error_rate.m for temporally interfering paths at receiver
% ber_info_no_interfere
%                 Structured array made and defined by ber_error_rate.m for direct paths only at receiver (no
%                 boundary-reflected path).
%
% ierr            0: no errors detected
%                 1: otherwise
%
%
% sim_c_info.mat  Matlab file with variables:   invars ber_info_interfere ber_info_no_interfere
%                 c3d_symbol_0 c3d_symbol_1  c3d_symbol_0_waveform c3d_symbol_1_waveform
%                 as defined in this program
% ber.fig .jpg    Figures showing bit error rates.
%
function [ber_info_interfere,ber_info_no_interfere,ierr]=simulate_c_information_measurement(invars)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preliminary Steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
progname='simulate_c_information_measurement.m';
% fprintf(1,'%s: See jlsb0.mat to debug. pausing\n',progname);save jlsb0;pause


dotest=1;
if(dotest==1)
    fprintf(1,'%s: Test case. Press ENTER to continue\n',progname);pause
    
    
    %-------------------------------------------------------------------------------------------------------------------
    % Parameters for simulate_two_waveforms
    %-------------------------------------------------------------------------------------------------------------------
    invars.f=3.948340440094052e+03;
    invars.tau=0.014;
    %invars.tau_after=[0.0001 0.0003 .0005];
    invars.tau_after=[0.02 0.06 0.1];
    invars.amp_big=1.5;
    invars.f_samp=1.d5;
    invars.n_tau_to_zero=3;
    invars.fign_ideal_symb=1;
    %-------------------------------------------------------------------------------------------------------------------
    % Paramters for bandpass filter
    %-------------------------------------------------------------------------------------------------------------------
    invars.stopband_freq1=3800;
    invars.passband_freq1=3900;
    invars.passband_freq2=4100;
    invars.stopband_freq2=4200;
    invars.stopband_atten1=30;
    invars.passband_ripple=1;
    invars.stopband_atten2=30;
    invars.design_method='ellip';
    invars.fign_filter=2;
    %-------------------------------------------------------------------------------------------------------------------
    % Paramters for bandpass filter at receiver
    %-------------------------------------------------------------------------------------------------------------------
    invars.receiver_filter.stopband_freq1=3800;
    invars.receiver_filter.passband_freq1=3900;
    invars.receiver_filter.passband_freq2=4100;
    invars.receiver_filter.stopband_freq2=4200;
    invars.receiver_filter.stopband_atten1=30;
    invars.receiver_filter.passband_ripple=1;
    invars.receiver_filter.stopband_atten2=30;
    invars.receiver_filter.design_method='ellip';
    invars.receiver_filter.fign_filter=2;
    %-------------------------------------------------------------------------------------------------------------------
    % Other parameters
    %-------------------------------------------------------------------------------------------------------------------
    
    invars.boundary_type=1;
    invars.c=1500;
    
    invars.method_interp='linear';
    invars.energy_loss=1;
    
    invars.snr=50;
    invars.n_symbol_transmit=400;
    invars.n_samps_before_c=10;
    invars.seedd=1;
    %-------------------------------------------------------------------------------------------------------------------
    % Parameters for bit error rate
    %-------------------------------------------------------------------------------------------------------------------


    invars.ber_thresh=0.1;

    invars.src_xz=[0 -11.746378838836542];
    invars.rec_xz=[12.697968050537270 -41.995196258872447];
end




% Initialize outputs
ierr=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute everyting wanted except bit error rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize random number generators
rng(invars.seedd);
% Do not override plot making in simulate_c_information_measurement_part1.m by setting make_plts to TRUE
make_plts=true(1,1);
verbose=1;

[h,ierr]=simulate_c_information_measurement_part1(invars,make_plts,verbose);
if(ierr~=0),return,end

if(h.t_direct > h.t_ideal_separate)
    fprintf(1,'Could distinguish between symbols 0 and 1 if infinite SNR and bandwidth.\n');
    fprintf(1,'h.t_ideal_separate is with respect to time optical  switch is thrown.\n');
    fprintf(1,'Pausing. Press ENTER to continue\n');
    pause
end

if(h.t_direct > h.t_ideal_separate)
    fprintf(1,'Could distinguish between symbols 0 and 1 if infinite SNR and bandwidth.\n');
    fprintf(1,'h.t_ideal_separate is with respect to time optical  switch is thrown.\n');
    fprintf(1,'Pausing. Press ENTER to continue\n');
    pause
end
if(~h.possible_faster_than_c)
    fprintf(1,'Press ENTER to compute BER anyway\n');
    pause
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derive bit-error rate (BER) for symbols at receivers for both direct and temporally inteferring paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------------------------------------------------------------------------------------------------------
% Derive bit error rate for temporally interferring signals at receiver: ber_info_interfere
%-----------------------------------------------------------------------------------------------------------------------
fign_ber=[];
[ber_info_interfere,ierr]=bit_error_rate(h.r0,h.r1,h.dt,h.t_s,h.ind_direct_na_arrive,...
    h.ind_tau_alpha,h.ind_na_arrive_no_filters_no_interference,h.snr,h.n_symbol_transmit,h.n_samps_before_c,h.d_rec,...
    fign_ber);
if(ierr~=0),return,end

%-----------------------------------------------------------------------------------------------------------------------
% Derive bit error rate for direct path only signals at receiver: ber_info_no_interfere
% This simulates scenario where only direct path is present
%-----------------------------------------------------------------------------------------------------------------------
fign_ber=[];
[ber_info_no_interfere,ierr]=bit_error_rate(h.r0_direct,h.r1_direct,h.dt,h.t_s,h.ind_direct_na_arrive,...
    h.ind_tau_alpha,h.ind_na_arrive_no_filters_no_interference,h.snr,h.n_symbol_transmit,h.n_samps_before_c,h.d_rec,...
    fign_ber);
if(ierr~=0),return,end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display BER results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(14);
clf;
semilogy(ber_info_interfere.tau,ber_info_interfere.ber,'k--');
hold on
semilogy(ber_info_no_interfere.tau,ber_info_no_interfere.ber,'k-');
xlabel('\tau (s)');
ylabel('BER');
title('solid: direct path dashed: interferring path \color{red} direct path no filters');


V=axis;
% Set axis limits
V(3)=0.0001;
V(4)=1;
axis(V);

% Plot time point of non-analyticity along direct path arrives in absence of any delays through filters and
% no temporal interference. It occurs at tau=0 according to bit_error_rate.m
% ber_info_interfere.tau=0 corresponds to ind_tau_alpha according to bit_error_rate.m
plot(repmat(h.t(h.ind_na_arrive_no_filters_no_interference)-h.t(h.ind_tau_alpha),1,2),V(3:4),'r-');

% Plot threshold
plot(V(1:2),repmat(invars.ber_thresh,1,2),'k-');
grid on



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save sim_c_info invars ber_info_interfere ber_info_no_interfere h

% Save BER figure
saveas(14,'ber','fig');
print -djpeg -r400 ber.jpg
