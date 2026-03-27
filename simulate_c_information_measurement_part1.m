% Is most of software from simulate_c_information_measurement.m before computing the bit error rate
%
% INPUTS
%
% invars          Structured array defined by simulate_c_information_measurement.m
% make_plts       1 x 1.  TRUE: render plots to screen
%                         FALSE: otherwise.  This overrides values in invars
% verbose         1 x 1.  0: Do not print out information to screen
%                         1: Otherwise
%
% OUTPUTS
%
% h               Structured array with fields
%
%  -r0            1 x nt.  r0(i) is signal at receiver due to symbol 0, before receiver's bandpsss, due to direct
%                 and reflected paths.
%                 r0 incluces effects due to finite bandwidth of optical switch.  r0(i) occurs at time t(i).
%  -r1            1 x nt.  Same as r0 except is due to symbol 1.
%  -r0_direct     1 x nt.  Same as r0 except due only to direct path.
%  -r1_direct     1 x nt.  Same as r1 "
%  -t             1 x nt.  Time axis (s) for above.  Defined by simulate_two_waveforms.m
%  -dt            1 x 1.   1/invars.f_samp.
%  -t_s           1 x 1.   Defined by bit_error_rate.m
%  -ind_direct_na_arrive
%                 1 x 1.  index of t corresponding to time soonest point of non-analyticity arrives at receiver for
%                 direct path.  The "soonest" means the soonest time derived from symbols 0 & 1 after each
%                 goes through the optical switch's filter.  The delay of each symbol through switch may differ.
%  -ind_tau_alpha 1 x 1.  Defined by bit_error_rate.m
%  -snr           1 x 1. Same as invars.snr.
%  -n_symbol_transmit
%                 1 x 1. Same as invars.n_symbol_transmit
%  -n_samps_before_c
%                 1 x 1. Same as invars.n_samps_before_c
%
%  -possible_faster_than_c
%                 1 x 1. Logical.  FALSE.  Impossible to transmit information about symbols faster than speed of
%                                          light (invars.c) in vacuum.
%                                  TRUE:   Might be possible to transmit information about symbols faster than
%                                          speed of light in vacuum.
%  -ind_na_arrive_no_filters_no_interference
%                  1 x 1. index of time axis, t, when point of non-analyticity arrives at receiver, excluding effects
%                  of optical switch and any bandpass effects of the receiver & assuming no temporal interference:
%  -ind_r0_r1_before_separate
%                  1 x 1.  largest index of r0 & r1 before they separate for first time
%  -crit_dist      1 x 1.  Critical distance (m) based on time scale of Gaussian envelopes
%  -crit_dist_optical_switch
%                  1 x 1.  Critical distance (m) based on bandwidth of optical switch
%  -c3d_symbol_0   1 x 1.  3D effective speed (m/s) for symbol 0  derived by
%                          1. cross-correlating symbol 0 at receiver (before receiver bandpass) with
%                             signal at source after optical switch
%                          2. Find lag of peak of absolute value of Hilbert transform of cross-correlation from 1.
%  -c3d_symbol_1   1 x 1.  Same as c3d_symbol_0 except uses symbol 1.
%
%  -c3d_symbol_0_waveform  3D effective speed (m/s) for symbol 0 derive by
%                          1. cross-correlating symbol 0 at receiver (before receiver bandpass) with
%                             signal at source after optical switch
%                          2. Deriving lag corresponding to max value of 1.
%
%  -c3d_symbol_1_waveform
%                  1 x 1.  Same as c3d_symbol_0_waveform except for symbol 1.
%  -t_ideal_separate
%                  1 x 1. Time (s) when symbols 0 and 1 first separate at receiver if filters have infinite bandwidth.
%                  This is elapsed time with respect to time when optical switch is thrown at source.
%  -t_direct       1 x 1. Time (s) of propagation of signal at speed invars.c between source and receiver. Equals
%                  distance between them divided by invars.c.
%  -d_rec          Structured array defined by output of design_bandpass_filter.m for filter parameters specified in
%                  invars.receiver_filter
%
% OUTPUTS
%
% ierr            0: no errors detected
%                 1: otherwise
%
function [h,ierr]=simulate_c_information_measurement_part1(invars,make_plts,verbose)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preliminary Steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
progname='simulate_c_information_measurement_part1.m';
% fprintf(1,'%s: See jlsb0.mat to debug. pausing\n',progname);save jlsb0;pause

% Initialize outputs
h=[];
ierr=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unpack inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------------------------------------------------
% Parameters for simulate_two_waveforms
%-------------------------------------------------------------------------------------------------------------------
f=invars.f;
tau=invars.tau;
tau_after=invars.tau_after;
amp_big=invars.amp_big;
f_samp=invars.f_samp;
n_tau_to_zero=invars.n_tau_to_zero;
fign_ideal_symb=invars.fign_ideal_symb;
%-------------------------------------------------------------------------------------------------------------------
% Paramters for bandpass filter for optical switch
%-------------------------------------------------------------------------------------------------------------------
stopband_freq1=invars.stopband_freq1;
passband_freq1=invars.passband_freq1;
passband_freq2=invars.passband_freq2;
stopband_freq2=invars.stopband_freq2;
stopband_atten1=invars.stopband_atten1;
passband_ripple=invars.passband_ripple;
stopband_atten2=invars.stopband_atten2;
design_method=invars.design_method;
fign_filter=invars.fign_filter;
%-------------------------------------------------------------------------------------------------------------------
% Paramters for bandpass filter at receiver
%-------------------------------------------------------------------------------------------------------------------
receiver_filter=invars.receiver_filter;
%-------------------------------------------------------------------------------------------------------------------
% Other parameters
%-------------------------------------------------------------------------------------------------------------------

boundary_type=invars.boundary_type;
c=invars.c;

method_interp=invars.method_interp;
energy_loss=invars.energy_loss;

snr=invars.snr;
n_symbol_transmit=invars.n_symbol_transmit;
n_samps_before_c=invars.n_samps_before_c;
%-------------------------------------------------------------------------------------------------------------------
% Parameters for source & receiver
%-------------------------------------------------------------------------------------------------------------------
src_xz=invars.src_xz;
rec_xz=invars.rec_xz;

ber_thresh=invars.ber_thresh;

% Derive approximate pulse resolution of emitted signal
pulse_resol=min([invars.tau invars.tau_after]);
% Derive critical distance from boundary: crit_dist
crit_dist=invars.c*pulse_resol/2;

% Derive critical distance due to bandwidth of optical switch
crit_dist_optical_switch=c/(passband_freq2-passband_freq1)/2;

if(verbose==1)
    fprintf(1,'%s: Critical distance based on time scales of Gaussian envelopes = %f m\n',progname,crit_dist);
    
    fprintf(1,' Critical distance based on bandwidth of optical switch = %f m\n',crit_dist_optical_switch);
    
    fprintf(1,'src_xz = (%20.10f, %20.10f) m\n',src_xz);
    fprintf(1,'rec_xz = (%20.10f, %20.10f) m\n',rec_xz);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create ideal emitted symbols
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(~make_plts)
    fign_ideal_symb=[];
end
[w,e,t,ierr]=simulate_two_waveforms(f,tau,tau_after,amp_big,f_samp,n_tau_to_zero,fign_ideal_symb);
if(ierr~=0),return,end

% Get # time samples each waveform
nt=length(t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Design bandpass filter to mimic emitted signal limited by analogy to optical switch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(~make_plts)
    fign_filter=[];
end
[d]=design_bandpass_filter(stopband_freq1,passband_freq1,passband_freq2,stopband_freq2,...
    stopband_atten1,passband_ripple,stopband_atten2,design_method,f_samp,fign_filter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply bandpass to ideal emitted signals to mimic affect of optical switch: y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[y]=filt_the_data(d,w,[nt nt]);

if(make_plts)
    % Plot ideal and bandpassed data
    figure(3);
    clf;
    subplot(2,1,1);
    plot(t,w(1,:),'k-',t,y(1,:),'r-');
    xlabel('s');
    ylabel('Symbol 0');
    title('Ideal \color{red} bandpassed \color{black} signal');
    grid on
    subplot(2,1,2);
    plot(t,w(2,:),'k-',t,y(2,:),'r-');
    xlabel('s');
    ylabel('Symbol 1');
    title('Ideal \color{red} bandpassed \color{black} signal');
    grid on
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derive delay through optical switch bandpass filter for symbols 0 & a: optical_delay_0 & optical_delay_1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find times >=0:
[loc]=find(t>=0);

% Cross correlate ideal emitted symbol 0 with same after optical delay
[z0_optical,lags_optical]=xcorr(y(1,loc),w(1,loc),'biased');
% Find lag corresponding to peak
[~,ind_peak_z0]=max(z0_optical);
% Get lag axis (s)
lags_optical_s=lags_optical/f_samp;

% Cross correlate ideal emitted symbol 1 with same after optical delay
[z1_optical,~]=xcorr(y(2,loc),w(2,loc),'biased');
% Find lag corresponding to peak
[~,ind_peak_z1]=max(z1_optical);

%-----------------------------------------------------------------------------------------------------------------------
% Derive delays due to optical switch on ideal emitted symbols 0 & 1
%-----------------------------------------------------------------------------------------------------------------------
optical_delay_0=lags_optical_s(ind_peak_z0);
optical_delay_1=lags_optical_s(ind_peak_z1);

%-----------------------------------------------------------------------------------------------------------------------
% Plot cross correlation between ideal and bandpassed symbols
%-----------------------------------------------------------------------------------------------------------------------
if(make_plts==1)
    figure(20);
    clf;
    % Plot z0
    subplot(2,1,1);
    plot(lags_optical_s,z0_optical,'k-',lags_optical_s(ind_peak_z0),z0_optical(ind_peak_z0),'r.');
    hold on;
    grid on
    V=axis;
    % Plot delay. loc is length n_loc
    plot(repmat(optical_delay_0,1,2),V(3:4),'g-');
    xlabel('Lag (s)');
    title(['x-corr ideal emitted sym 0 with same after optical switch. \color{green} delay= ' num2str(optical_delay_0)]);
    
    % Plot z1
    subplot(2,1,2);
    plot(lags_optical_s,z1_optical,'k-',lags_optical_s(ind_peak_z1),z1_optical(ind_peak_z1),'r.');
    hold on;
    grid on
    V=axis;
    % Plot delay
    plot(repmat(optical_delay_1,1,2),V(3:4),'g-');
    xlabel('Lag (s)');
    title(['x-corr ideal emitted sym 1 with same after optical switch. \color{green} delay= ' num2str(optical_delay_1)]);
end % if(make_plts==1)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derive distances and propagation times of direct and boundary reflected paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derive distance of direct path
d_direct=sqrt( (src_xz(1) - rec_xz(1))^2 + (src_xz(2)-rec_xz(2))^2);
% Derive delay along direct path
t_direct=d_direct/c;


% Derive distance of reflected path from soruce to receiver: d_refl
d_refl=sqrt( (src_xz(1) - rec_xz(1))^2 + (-src_xz(2)-rec_xz(2))^2);
% Derive delay along reflected path
t_refl=d_refl/c;

%-----------------------------------------------------------------------------------------------------------------------
% Derive approximate pulse resolution of emitted signal
%-----------------------------------------------------------------------------------------------------------------------
pulse_resol=min([tau tau_after]);

if(verbose==1)
    % Tell user if source is within cricitical distance from boundary
    if(pulse_resol< t_refl - t_direct)
        fprintf(1,'%s: pulse_resol = %f s < difference in time between reflected and direct paths = %f s\n',progname,...
            pulse_resol,t_refl-t_direct);
    else
        fprintf(1,'%s: pulse resol = %f s. t_refl-t_direct = %f s.\n',progname,...
            pulse_resol,t_refl-t_direct);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derive time axis for simulation: t_axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% derive time last part of reflected signal arrives at receiver: t_last
t_last=max(t)+t_refl;
% Derive number samples after max(t) to get to t_last: n
n=ceil((t_last-max(t))*f_samp)+1;
% Derive time axis values after t_max
t_after=max(t) + (1:n)/f_samp;
% Make taxis for experiment
t_axis=[t t_after];

% Get length of t_axis
n_t_axis=length(t_axis);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define y for all times on t_axis: yy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% yy is bandpassed signal from source where bandpass mimics optical switch
yy=zeros(2,n_t_axis);
yy(:,1:nt)=y;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate both signals (symbols) at receiver for both signals before including effects of receiver bandwidth:
% r0_before_rec_filt and r1_before_rec_filt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------------------------------------------------------------------------------------------------------
% Get symbol 0 at receiver
%-----------------------------------------------------------------------------------------------------------------------
% symbol 0 along direct path: r0_direct
if(energy_loss==1)
    r0_direct=interp1(t_axis,yy(1,:),t_axis-t_direct,method_interp,0)/d_direct;
else
    r0_direct=interp1(t_axis,yy(1,:),t_axis-t_direct,method_interp,0);
end

% symbol 0 along reflected path: r0_reflect
if(energy_loss==1)
    r0_reflect=interp1(t_axis,yy(1,:),t_axis-t_refl,method_interp,0)/d_refl;
else
    r0_reflect=interp1(t_axis,yy(1,:),t_axis-t_refl,method_interp,0);
end
if(boundary_type==1)
    % No pressure fluctuation at surface.  Phase flips by 180 degreees
    r0_reflect=-r0_reflect;
else
    % Nothing to modify
end

% Plot emitted and direct path if wanted
if(1==2)
    figure(3);
    clf;
    ax1=subplot(2,1,1);
    plot(ax1,t_axis,yy(1,:),'k-');
    title(ax1,'Emitted symbol 0');
    ax2=subplot(2,1,2);
    plot(ax2,t_axis,r0_direct,'k-');
    title(ax2,'Direct Path');
    grid on
    % Link axes
    linkaxes([ax1 ax2],'xy');
end
%-----------------------------------------------------------------------------------------------------------------------
% Get symbol 1 at receiver
%-----------------------------------------------------------------------------------------------------------------------
% symbol 1 along direct path: r1_direct
if(energy_loss==1)
    r1_direct=interp1(t_axis,yy(2,:),t_axis-t_direct,method_interp,0)/d_direct;
else
    r1_direct=interp1(t_axis,yy(2,:),t_axis-t_direct,method_interp,0);
end

% symbol 1 along reflected path: r1_reflect
if(energy_loss==1)
    r1_reflect=interp1(t_axis,yy(2,:),t_axis-t_refl,method_interp,0)/d_refl;
else
    r1_reflect=interp1(t_axis,yy(2,:),t_axis-t_refl,method_interp,0);
end

if(boundary_type==1)
    % No pressure fluctuation at surface.  Phase flips by 180 degreees
    r1_reflect=-r1_reflect;
else
    % Nothing to modify
end

% Plot emitted and direct path if wanted
if(1==2)
    figure(4);
    clf;
    ax1=subplot(2,1,1);
    plot(ax1,t_axis,yy(2,:),'k-');
    title(ax1,'Emitted symbol 1');
    ax2=subplot(2,1,2);
    plot(ax2,t_axis,r1_direct,'k-');
    title(ax2,'Direct Path');
    grid on
    % Link axes
    linkaxes([ax1 ax2],'xy');
end

%-----------------------------------------------------------------------------------------------------------------------
% Derive signals at receiver before bandpassing at receiver: r0 and r1
%-----------------------------------------------------------------------------------------------------------------------
r0=r0_direct + r0_reflect;
r1=r1_direct + r1_reflect;

%-----------------------------------------------------------------------------------------------------------------------
% Show symbol 0 decomposition before bandpassing at receiver
%-----------------------------------------------------------------------------------------------------------------------
if(make_plts)
    figure(5);
    clf;
    ax1=subplot(3,1,1);
    plot(ax1,t_axis,r0_direct);
    title(ax1,'Received symbol 0 direct path before rec. bandpass');
    xlabel('s');
    grid on
    
    ax2=subplot(3,1,2);
    plot(ax2,t_axis,r0_reflect);
    title(ax2,'Received symbol 0 refl path before rec. bandpass');
    xlabel('s')
    grid on
    
    ax3=subplot(3,1,3);
    plot(ax3,t_axis,r0);
    title(ax3,'Received symbol 0 before rec. bandpass');
    xlabel('s');
    % Link axes
    linkaxes([ax1 ax2 ax3 ],'xy');
    
    
    %-----------------------------------------------------------------------------------------------------------------------
    % Show symbol 1 decomposition before bandpassing at receiver
    %-----------------------------------------------------------------------------------------------------------------------
    figure(6);
    clf;
    ax1=subplot(3,1,1);
    plot(ax1,t_axis,r1_direct);
    title(ax1,'Received symbol 1 direct path before rec. bandpass');
    xlabel('s');
    grid on
    
    ax2=subplot(3,1,2);
    plot(ax2,t_axis,r1_reflect);
    title(ax2,'Received symbol 1 refl path before rec. bandpass');
    xlabel('s')
    grid on
    ax3=subplot(3,1,3);
    plot(ax3,t_axis,r1);
    title(ax3,'Received symbol 1 before rec. bandpass');
    xlabel('s')
    grid on
    % Link axes
    linkaxes([ax1 ax2 ax3 ],'xy');
    
    %-----------------------------------------------------------------------------------------------------------------------
    % Show symbols 0 & 1  before bandpassing at receiver
    %-----------------------------------------------------------------------------------------------------------------------
    figure(7);
    clf;
    plot(t_axis,r0,'k-',t_axis,r1,'g-');
    title('Symbols 0 and \color{green} 1 \color{black} at receiver before bandpass');
    xlabel('s');
    grid on
    
end % if(make_plts)
% Save signals arriving at receiver before receiver's bandpass: r0_before_rec_filt & r1_before_rec_filt
r0_before_rec_filt=r0;
r1_before_rec_filt=r1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Design bandpass filter to mimic receiver's finite bandwidth: d_rec
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(~make_plts)
    receiver_filter.fign_filter=[];
end
[d_rec]=design_bandpass_filter(receiver_filter.stopband_freq1,receiver_filter.passband_freq1,...
    receiver_filter.passband_freq2,receiver_filter.stopband_freq2,...
    receiver_filter.stopband_atten1,receiver_filter.passband_ripple,receiver_filter.stopband_atten2,...
    receiver_filter.design_method,f_samp,receiver_filter.fign_filter);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply bandpass to signals arriving at receiver: overwrite r0 & r1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[r_filt]=filt_the_data(d_rec,[r0_before_rec_filt;r1_before_rec_filt],[n_t_axis n_t_axis]);
% Overwrite r0 and r1
r0=r_filt(1,:);
r1=r_filt(2,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derive delay through receiver's bandpass filter for symbols 0 & 1: receiver_delay_0 & receiver_delay_1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cross correlate ideal emitted symbol 0 with same after optical delay
[z0_rec,lags_rec]=xcorr(r0(loc),r0_before_rec_filt(loc),'biased');

% Find lag corresponding to peak
[~,ind_peak_z0_rec]=max(z0_rec);
% Get lag axis (s)
lags_rec_s=lags_rec/f_samp;

% Cross correlate ideal emitted symbol 1 with same after optical delay
[z1_rec,~]=xcorr(r1(loc),r1_before_rec_filt(loc),'biased');
% Find lag corresponding to peak
[~,ind_peak_z1_rec]=max(z1_rec);

%-----------------------------------------------------------------------------------------------------------------------
% Derive delays due to optical switch on ideal emitted symbols 0 & 1
%-----------------------------------------------------------------------------------------------------------------------
receiver_delay_0=lags_rec_s(ind_peak_z0_rec);
receiver_delay_1=lags_rec_s(ind_peak_z1_rec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get absolute value Hilbert transform of both symbols at receiver after accounting for bandwidth of receiver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ah_r0=abs(hilbert(r0));
ah_r1=abs(hilbert(r1));

if(make_plts)
    figure(9);
    clf;
    plot(t_axis,ah_r0,'k-',t_axis,ah_r1,'g-');
    title('Abs Hilbert Symbol 0 and \color{green} 1 \color{black} at receiver after its bandpass');
    xlabel('s');
    grid on
end % if(make_plts)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derive c3d by cross-correlating y(1,:) with r0 and y(2,:) with r1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add trailing zeros to y so same length as r0 and r1: y_pad0 & y_pad1
[y_pad0,~]=make_same_length(y(1,:),r0);
[y_pad1,~]=make_same_length(y(2,:),r1);

% cross correlate received signal for symbol 0 with emitted symbol 0
[z0,lag_axis]=xcorr(r0,y_pad0,'biased');
% cross correlate received signal for symbol 1 with emitted symbol 1
[z1,~]       =xcorr(r1,y_pad1,'biased');

% Convert lag axis to seconds
lag_axis=lag_axis/f_samp;

% Derive delay of waveform, r0, with respect to  waveform, y_pad0: pk_lag_0_waveform.
% Take abs value of z0 because lag whee abs(z0) is largest may have negative value.
[~,ind_pk_0_waveform]=max(abs(z0));
pk_lag_0_waveform=lag_axis(ind_pk_0_waveform);
% Derive delay of waveform, r1, with respect to  waveform, y_pad1: pk_lag_1_waveform
% Take abs value of z1 because lag whee abs(z1) is largest may have negative value.
[~,ind_pk_1_waveform]=max(abs(z1));
pk_lag_1_waveform=lag_axis(ind_pk_1_waveform);

% Get abs value of Hilbert transforms
ah_z0=abs(hilbert(z0));
ah_z1=abs(hilbert(z1));

% Find peaks
[~,ind_pk_0]=max(ah_z0);
[~,ind_pk_1]=max(ah_z1);

% Get lags (s). Positive means cross correlation peak later at receiver than source
pk_lag_0=lag_axis(ind_pk_0);
pk_lag_1=lag_axis(ind_pk_1);


%-----------------------------------------------------------------------------------------------------------------------
% Plot abs value Hilbert transforms near peak
%-----------------------------------------------------------------------------------------------------------------------
%  Set limits for plots
% abs_tau_max_plt=0.0001;
% V=[-abs_tau_max_plt abs_tau_max_plt 0 2];

if(make_plts)
    % Plot abs value Hilbert transform ccf and lag of peak
    % Set linewidth
    lw=2;
    % Allocate memory for axis limits
    V=zeros(2,4);
    figure(15);
    clf;
    subplot(2,1,1);
    plot(lag_axis,ah_z0,'k-',lag_axis(ind_pk_0),ah_z0(ind_pk_0),'r.','markersize',7,'linewidth',lw);
    %axis(V);
    % Store axis limits
    V(1,:)=axis;

    hold on;
    % Plot time of direct arrival
    plot(repmat(t_direct,1,2),V(1,3:4),'k--');
    grid on
    xlabel('Lag (s)');
    title('Abs Hilbert CCF Symbol 0. \color{red} peak \color{black} Lag Direct Path');
    subplot(2,1,2);
    plot(lag_axis,ah_z1,'k-',lag_axis(ind_pk_1),ah_z1(ind_pk_1),'r.','markersize',7,'linewidth',lw);
    % Store axis limits
    V(2,:)=axis;
    hold on;

    grid on;
    xlabel('Lag (s)');
    title('Abs Hilbert CCF Symbol 1. \color{red} peak \color{black} Lag Direct Path');

    % Get limits to use for both subplots: V_use
    V_use=[min(V(:,1)) max(V(:,2)) min(V(:,3)) max(V(:,4))];
    % Make axis limits same and draw time direct path
    for i=1:2
        subplot(2,1,i);
        axis(V_use);
        % Plot time of direct arrival
        plot(repmat(t_direct,1,2),V_use(3:4),'k--');
    end
    
end % if(make_plts)

%-----------------------------------------------------------------------------------------------------------------------
% Derive c3d for both symbols
%-----------------------------------------------------------------------------------------------------------------------
% Derive c3d if using peak of abs value of Hilbert transform
c3d_symbol_0=d_direct/pk_lag_0;
c3d_symbol_1=d_direct/pk_lag_1;

% Derive c3d if using peak of cross-correlated waveforms
c3d_symbol_0_waveform=d_direct/pk_lag_0_waveform;
c3d_symbol_1_waveform=d_direct/pk_lag_1_waveform;

if(verbose==1)
    fprintf(1,'\n');
    fprintf(1,'c3d for symbol 0 and 1 are %f and %f m/s respectively based on peak of abs Hilbert\n',...
        c3d_symbol_0,c3d_symbol_1);
    fprintf(1,'c3d for symbol 0 and 1 are %f and %f m/s respectively based on peak of cross-correlation waveforms\n',...
        c3d_symbol_0_waveform,c3d_symbol_1_waveform);
    fprintf(1,'\n');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show symbols 0 & 1 compared with only their direct signals after bandpassing at receiver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(make_plts)
    %-------------------------------------------------------------------------------------------------------------------
    % Plot symbol r0 compared with its direct path at receiver before bandpassing at receiver
    %-------------------------------------------------------------------------------------------------------------------
    % Set linewidth
    lw=2;
    figure(10)
    clf;
   
    % Show waveforms
    ax1=subplot(2,1,1);
    plot(ax1,t_axis,r0,'g-',t_axis,r0_direct,'k-','linewidth',lw);
    title(ax1,'r0 direct \color{green} r0  at rec after bandpass');
    xlabel(ax1,'s');
    grid on
    
    % Show abs Hilbert transform waveforms
    ax2=subplot(2,1,2);
    plot(ax2,t_axis,ah_r0,'g-',t_axis,abs(hilbert(r0_direct)),'k-','linewidth',lw);
    title(ax2,'Abs Hilbert transform. Dashed: Time ideal non-analytic direct path');
    xlabel(ax2,'s');
    % Link axes
    linkaxes([ax1 ax2],'xy');
    grid on
    %-------------------------------------------------------------------------------------------------------------------
    % Plot symbol r1 compared with its direct path at receiver before bandpassing at receiver
    %-------------------------------------------------------------------------------------------------------------------
    figure(11)
    clf;
    % Show waveforms
    ax1=subplot(2,1,1);
    plot(ax1,t_axis,r1,'k-',t_axis,r1_direct,'g-','linewidth',lw);
   
    title(ax1,'r1 direct \color{green} r1  at rec after bandpass');
    xlabel(ax1,'s');
    grid on
    
    % Show abs Hilbert transform waveforms
    ax2=subplot(2,1,2);
    plot(ax2,t_axis,ah_r1,'g-',t_axis,abs(hilbert(r1_direct)),'k-','linewidth',lw);
    title(ax2,'Abs Hilbert transform. Dashed: Time ideal non-analytic direct path');
    xlabel(ax2,'s');
    % Link axes
    linkaxes([ax1 ax2],'xy');
    grid on;
end % if(make_plts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find point in time just before r0 & r1 start to separate: ind_r0_r1_before_separate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get difference in temporally arriving signals: dr
dr=r1-r0;

% Find first index of abs(dr)>0
[ind_separate]=find(abs(dr)>0,1);
% Get largest index of r0 & r1 before they separate for first time
ind_r0_r1_before_separate=ind_separate-1;

if(make_plts)
    figure(12);
    clf
    plot(t_axis,r0,'k-',t_axis,r1,'g-');
    V=axis;
    % Plot point of separation
    hold on;
    plot(repmat(t_axis(ind_r0_r1_before_separate),1,2),V(3:4),'r-');
    xlabel('(s)');
    title('symbol 0 and \color{green} 1 \color{black} at receiver due to direct+reflected paths. \color{red} separate');
end % if(make_plts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply bandpass to only direct-path signals arriving at receiver: Overwrite r0_direct & r1_direct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[r_direct_filt]=filt_the_data(d_rec,[r0_direct;r1_direct],[n_t_axis n_t_axis]);
% Overwrite r0 and r1
r0_direct=r_direct_filt(1,:);
r1_direct=r_direct_filt(2,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find point in time just before r0_direct & r1_direct start to separate: ind_r0_r1_direct_before_separate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get difference in temporally arriving signals: dr
dr=r1_direct-r0_direct;

% Find first index of abs(dr)>0
[ind_separate]=find(abs(dr)>0,1);
% Get largest index of r0 & r1 before they separate for first time
ind_r0_r1_direct_before_separate=ind_separate-1;

if(make_plts)
    figure(13);
    clf
    plot(t_axis,r0_direct,'k-',t_axis,r1_direct,'g-');
    V=axis;
    % Plot point of separation
    hold on;
    plot(repmat(t_axis(ind_r0_r1_direct_before_separate),1,2),V(3:4),'r-');
    xlabel('(s)');
    title('direct paths symbol 0 and \color{green} 1 \color{black} at receiver due to direct paths. \color{red} separate');
end % if(make_plts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derive information about symbols at receivers for both direct and temporally inteferring paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BER for direct path mimics scenario in vacuum
% "       symbol 0 (r0) mimics scenario of fast light propagation
%-----------------------------------------------------------------------------------------------------------------------
% Prepare some inputs for bit_error_rate.m
%-----------------------------------------------------------------------------------------------------------------------
dt=1/f_samp;
t_s=0;
%-----------------------------------------------------------------------------------------------------------------------
% Derive index of time axis, t, when point of non-analyticity arrives at receiver, excluding effects
% of optical switch and any bandpass effects of the receiver & assuming no temporal interference:
% ind_na_arrive_no_filters_no_interference
%-----------------------------------------------------------------------------------------------------------------------
% t(ind_na_arrive_no_filters) is time the point of non-analyticity arrives at receiver in vacuum
% if there is only the direct path and there are not bandpass filters.    t(loc(1)) is time the
% switch is thrown to transmit symbol 0 or 1.  If transmitted signal is electromagnetic wave, then
% t(ind_na_arrive_no_filters) - t(loc(1)) is the time the point of non-analyticity arrives if propagating
% through a vacuum and there is no temporal interference.
ind_na_arrive_no_filters_no_interference=round(t_direct*f_samp) + loc(1);

%-----------------------------------------------------------------------------------------------------------------------
% Determine if it is possible the symbols can be detected at receiver sooner than the time it takes
% light to go from time of optical switch to time the point of non-analyticity arrives at receiver
%-----------------------------------------------------------------------------------------------------------------------
% If symbols 0 and 1 are identical at times equal to or greater than the time needed for information to travel
% and the speed of light, there is no possibility their identities cannot be detected at times
% less than that needed to transmit info at the speed of light (c).
if(ind_r0_r1_before_separate+1 < ind_na_arrive_no_filters_no_interference)
    % It may be possible to detect which symbol is transmitted faster than light speed
    possible_faster_than_c=true(1,1);
else
    possible_faster_than_c=false(1,1);
end

if(verbose==1)
    if(possible_faster_than_c)
        % It is possible the detection of symbol 0 and 1 can be made faster than the speed of light
        fprintf(1,'%s: It might be possible to distinguish symbol 0 from 1 at times less than\n',progname);
        fprintf(1,' information traveling at the speed of light (c) in vacuum.\n');
    else
        fprintf(1,'%s: Not possible to distinguish between symbols 0 and 1 at times less than \n',progname);
        fprintf(1,' information traveling at the speed of light (c)\n');
        
    end
end
% Derive the least possible time needed to distinguish symbols 0 & 1:
t_to_distinguish=t_axis(ind_r0_r1_before_separate+1);
% Derive time needed to transmit information at speed of light (c) in vacuum: t_distinguish_light_speed
t_distinguish_light_speed=t_axis(ind_na_arrive_no_filters_no_interference);

%-----------------------------------------------------------------------------------------------------------------------
% Derive ind_direct_na_arrive
%-----------------------------------------------------------------------------------------------------------------------
% Derive time when point of non-analyticity arrives at receiver for direct path (after receiver bandwidth if modeled)
% for symbols 0 & 1
t0_arrives_na=t_direct+optical_delay_0 + receiver_delay_0;
t1_arrives_na=t_direct+optical_delay_1 + receiver_delay_1;

% Compute index of r0_direct when its point of non-analyticity arrives at receiver
ind0=ceil(t0_arrives_na*f_samp);
ind1=ceil(t1_arrives_na*f_samp);

% loc(1) corresponds to t=0 for e, i.e. e(:,loc) corresponds to t=0.  This is time the optical switch is thrown.
ind_direct_na_arrive=loc(1) + min([ind0 ind1]);
% Adjust ind_direct_na_arrive if needed till r0_direct(ind_direct_na_arrive) = r1_direct(ind_direct_na_arrive)
i_same=NaN;
for i=ind_direct_na_arrive:-1:1
    if(r0_direct(i)==r1_direct(i))
        % r0_direct(i) = r1_direct(i)
        i_same=i;
        break
    end
end

if(isnan(i_same))
    fprintf(1,'%s: BUG. Could not find and index, i, where r0_direct(i)-r1_direct(i)\n',progname);
    ierr=1;
    return
else
    % r0_direct(ind_direct_na_arrive) = r1_direct(ind_direct_na_arrive)
    ind_direct_na_arrive=i_same;
end

if(make_plts)
    %-------------------------------------------------------------------------------------------------------------------
    % Plot time point of non-analyticity arrives for direct path on figs 10 & 11. This is time if no filters anywhere
    %-------------------------------------------------------------------------------------------------------------------
    for i=10:11
        figure(i);
        subplot(2,1,2);
        V=axis;

        hold on
        plot(repmat(t_direct,1,2),V(3:4),'k--');
    end

    %-------------------------------------------------------------------------------------------------------------------
    % Get index to use for bit_error_rate indicating point in time before symbols separate: ind_tau_alpha
    %-------------------------------------------------------------------------------------------------------------------
    % Want to use same value of ind_tau_alpha for both calls to bit_error_rate.m so bit error rates have the same
    % starting values of tau
    figure(9);
    V=axis;
    hold on
    plot(repmat(t_direct,1,2),V(3:4),'k--');
end % if(make_plts)


%-----------------------------------------------------------------------------------------------------------------------
% Plot time point of non-analyticity arrives for direct path on figs 9. This is time if no filters anywhere
%-----------------------------------------------------------------------------------------------------------------------
ind_tau_alpha=min([ind_r0_r1_before_separate ind_r0_r1_direct_before_separate]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derive time when ideal symbols separate at receiver if all  filters have infinite bandwidth: t_ideal_separate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------------------------------------------------------------------------------------------------------
% Get ideal symbol 0 at receiver: r0_ideal
%-----------------------------------------------------------------------------------------------------------------------
% symbol 0 along direct path: r0_direct_ideal
if(energy_loss==1)
    r0_direct_ideal=interp1(t,e(1,:),t-t_direct,method_interp,0)/d_direct;
else
    r0_direct_ideal=interp1(t_e,e(1,:),t-t_direct,method_interp,0);
end

% symbol 0 along reflected path: r0_reflect_ideal
if(energy_loss==1)
    r0_reflect_ideal=interp1(t,e(1,:),t-t_refl,method_interp,0)/d_refl;
else
    r0_reflect_ideal=interp1(t,e(1,:),t-t_refl,method_interp,0);
end
if(boundary_type==1)
    % No pressure fluctuation at surface.  Phase flips by 180 degreees
    r0_reflect_ideal=-r0_reflect_ideal;
else
    % Nothing to modify
end

r0_ideal=r0_direct_ideal + r0_reflect_ideal;

%-----------------------------------------------------------------------------------------------------------------------
% Get ideal symbol 0 at receiver: r1_ideal
%-----------------------------------------------------------------------------------------------------------------------
% symbol 1 along direct path: r0_direct_ideal
if(energy_loss==1)
    r1_direct_ideal=interp1(t,e(2,:),t-t_direct,method_interp,0)/d_direct;
else
    r1_direct_ideal=interp1(t_e,e(2,:),t-t_direct,method_interp,0);
end

% symbol 1 along reflected path: r0_reflect_ideal
if(energy_loss==1)
    r1_reflect_ideal=interp1(t,e(2,:),t-t_refl,method_interp,0)/d_refl;
else
    r1_reflect_ideal=interp1(t,e(2,:),t-t_refl,method_interp,0);
end
if(boundary_type==1)
    % No pressure fluctuation at surface.  Phase flips by 180 degreees
    r1_reflect_ideal=-r1_reflect_ideal;
else
    % Nothing to modify
end

r1_ideal=r1_direct_ideal + r1_reflect_ideal;


%-----------------------------------------------------------------------------------------------------------------------
% Identify first index of r0_ideal and r1_ideal  when separation first occurs: ind_ideal_separate
%-----------------------------------------------------------------------------------------------------------------------
% t is 1-1 with r1_ideal andd r0_ideal
[ind_ideal_separate]=find(r1_ideal-r0_ideal~=0,1);
% Get associated time with respect to time switch is thrown
t_ideal_separate=t(ind_ideal_separate);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h.r0=r0;
h.r1=r1;
h.r0_direct=r0_direct;
h.r1_direct=r1_direct;
h.t=t_axis;
h.dt=dt;
h.t_s=t_s;
h.ind_direct_na_arrive=ind_direct_na_arrive;
h.ind_tau_alpha=ind_tau_alpha;
h.snr=snr;
h.n_symbol_transmit=n_symbol_transmit;
h.n_samps_before_c=n_samps_before_c;
h.possible_faster_than_c=possible_faster_than_c;
h.ind_na_arrive_no_filters_no_interference=ind_na_arrive_no_filters_no_interference;
h.ind_r0_r1_before_separate=ind_r0_r1_before_separate;
h.crit_dist=crit_dist;
h.crit_dist_optical_switch=crit_dist_optical_switch;
h.c3d_symbol_0=c3d_symbol_0;
h.c3d_symbol_1=c3d_symbol_1;
h.c3d_symbol_0_waveform=c3d_symbol_0_waveform;
h.c3d_symbol_1_waveform=c3d_symbol_1_waveform;
h.t_ideal_separate=t_ideal_separate;
h.t_direct=t_direct;
h.d_rec=d_rec;






