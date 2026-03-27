% Runs simulate_c_information_measurement_part1.m for many cases to see if any might be able to transmit
% information faster than speed of light in vacuum.
%
% INPUTS
%
% invars           Structured array defined by simulate_c_information_measurement.m
% make_plts        1 x 1.  TRUE: render plots to screen from simulate_c_information_measurement_part1.m
%                          FALSE: otherwise.  This overrides values in invars
%
% OUTPUTS
%   
% hh               Structured array with fields
% ierr             0: no errors detected
%                  1: otherwise
%
%
function [hh,ierr]=many_cases_simulate_c_information_measurement_part1(invars,make_plts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preliminary Steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
progname='many_cases_simulate_c_information_measurement_part1.m';
% fprintf(1,'%s: See jlsb0.mat to debug. pausing\n',progname);save jlsb0;pause



% Initialize outputs
dotest=1;
if(dotest==1)
    fprintf(1,'%s: Test case. Press ENTER to continue\n',progname);pause
    
    %-------------------------------------------------------------------------------------------------------------------
    % Parameters for simulate_two_waveforms
    %-------------------------------------------------------------------------------------------------------------------
    invars.f=3950;
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
    

end % if(dotest==1)



% Initialize outputs
hh=[];
ierr=0;


%-----------------------------------------------------------------------------------------------------------------------
% Parameters for many cases
%-----------------------------------------------------------------------------------------------------------------------
% Specify number cases to run
n_cases=100;
% Specify which critical distance to use for setting bounds of source and receiver: which_crit
% which_crit:  1:  use crit_dist
%              2:  use critdist_optical_switch
which_crit=1;

% Choose multiplication factor, f_mult, for setting  bounds of source & receiver location: f_mult
f_mult=4;

% Set bounds for center frequency of signal (Hz)
f_bnds=[3940 3960];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run simulate_c_information_measurement_part1.m for many cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Override plot making in simulate_c_information_measurement_part1.m by setting make_plts to TRUE
make_plts=false(1,1);
% Suppress messages to screen from simulate_c_information_measurement_part1.m
verbose=0;

% Derive spatial bounds for source and receiver: loc_bnds
[loc_bnds,crit_dist_use,ierr]=derive_loc_bnds(invars.tau,invars.tau_after,invars.passband_freq1,...
    invars.passband_freq2,...
    which_crit,f_mult,invars.c);
if(ierr~=0),return,end

% Allocate memory to store c3d for each case. All defined by simulate_c_information_measurement_part1.m
c3d_symbol_0=nan(1,n_cases);
c3d_symbol_1=nan(1,n_cases);
c3d_symbol_0_waveform=nan(1,n_cases);
c3d_symbol_1_waveform=nan(1,n_cases);
ind_na_arrive_no_filters_no_interference=nan(1,n_cases);
ind_r0_r1_before_separate=nan(1,n_cases);
% set interval to report progress to screen




report_interval=100;
t_start=tic;

for icase=1:n_cases
    if(rem(icase,report_interval)==0)
        % Get elapsed time so far
        t_elap=toc(t_start);
        % Get average time per case
        t_per_case=t_elap/icase;
        % Predict time till done
        t_till_done=(n_cases-icase)*t_per_case;
        fprintf(1,'Start case %d of %d. Done in %7.1f min.\n',icase,n_cases,t_till_done/60);
    end
    % Initialize random number generators for this case
    rng(invars.seedd + icase - 1);
    
    % Choose center frequency
    invars.f=f_bnds(1)+rand(1,1)*diff(f_bnds);
    
    % Choose z coordinate of source
    src_z=rand(1,1)*diff(loc_bnds(2,:)) + loc_bnds(2,1);
    % Choose z coordinate of receiver
    rec_z=rand(1,1)*diff(loc_bnds(2,:)) + loc_bnds(2,1);
    % Choose x coordinate of receiver
    rec_x=rand(1,1)*diff(loc_bnds(1,:)) + loc_bnds(1,1);
    
    % Insert source/receiver locations in invars
    invars.src_xz(2)=src_z;
    invars.rec_xz=[rec_x rec_z];
    
    [h,ierr]=simulate_c_information_measurement_part1(invars,make_plts,verbose);
    if(ierr~=0),return,end
    
    % Store some outputs
    ind_na_arrive_no_filters_no_interference(icase)=h.ind_na_arrive_no_filters_no_interference;
    ind_r0_r1_before_separate(icase)=h.ind_r0_r1_before_separate;
    c3d_symbol_0(icase)=h.c3d_symbol_0;
    c3d_symbol_1(icase)=h.c3d_symbol_1;
    c3d_symbol_0_waveform(icase)=h.c3d_symbol_0_waveform;
    c3d_symbol_1_waveform(icase)=h.c3d_symbol_1_waveform;
    
    if(h.possible_faster_than_c)
       fprintf(1,'%s: Found case where information might be transmitted faster than c. Pausing\n',progname);
       pause
    end

    if(h.t_direct > h.t_ideal_separate)
        fprintf(1,'Could distinguish between symbols 0 and 1 if infinite SNR and bandwidth.\n');
        fprintf(1,'h.t_ideal_separate is with respect to time optical  switch is thrown.\n');
        fprintf(1,'Pausing. Press ENTER to continue\n');
        pause
    end

end % for icases=1:n_cases

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set # subplots
nsub=2;
% Allocate memory to store handle for each subplot
    g.ax=[];
    ax_all=repmat(g,1,nsub);
    clear g
figure(1);
clf;
%-----------------------------------------------------------------------------------------------------------------------
% Plot c3d speeds
%-----------------------------------------------------------------------------------------------------------------------
ax_all(1).ax=subplot(nsub,1,1);
plot(ax_all(1).ax,1:n_cases,c3d_symbol_0,'k.',1:n_cases,c3d_symbol_1,'g.');
ylabel(ax_all(1).ax,'c3d (m/s)');
title(ax_all(1).ax,'Symbol 0 hilbert \color{green} 1 hilbert');

%-----------------------------------------------------------------------------------------------------------------------
% Plot min time needed to transmit information compared to info transmitted at light speed (invars.c)
%-----------------------------------------------------------------------------------------------------------------------

% Get min possible time needed to transmit informaiton: t_min
t_min=h.t(ind_r0_r1_before_separate+1);
% Get time to transmit information at light speed
t_light=h.t(ind_na_arrive_no_filters_no_interference);
% Get difference
dt=t_min-t_light;
ax_all(2).ax=subplot(nsub,1,2);
plot(ax_all(2).ax,1:n_cases,dt,'k.');
hold on;

grid on
xlabel(ax_all(2).ax,'Case');
ylabel(ax_all(2).ax,'s');
title(ax_all(2).ax,'Min time to transmit info - light speed info time');


% Link axes
% linkaxes([ax1 ax2],'x');
s='linkaxes([';
for i=1:nsub
    s=[s 'ax_all(' num2str(i) ').ax '];
end
s=[s '],''x'');'];
eval(s);

