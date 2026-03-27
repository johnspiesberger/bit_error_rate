% Computes exact time-domain solution of the linear acoustic wave equation
% in water of constant speed of sound with direct path and path reflecting from flat boundary.
% Boundary is either 0 pressure (surface of ocean) or hard boundary (like rocks).
%
% INPUTS
%
% c             1 x 1.  speed of sound in water (m/s)
% src_xz        1 x 2.  x = src_xz(1) (m) and z=src_xz (m) are the horizontal and vertical coordinates of a source.
%                       z is 0 at ocean's surface and is positive up.  source also located at y=0,  in 3D space.
%                       z is 0 at ocean's surface and is positive up.  source also located at y=0,  in 3D space.
% rec_xz        1 x 2. Is the x and z Cartesian coordinates of receiver
% tau           1 x 1. emitted pulse is exp(-(t/tau)^2)*cos(omega*t) with tau in seconds and t in seconds
% omega         1 x 1. See tau above. Units are rad/s
% g             1 x 1. Sample frequency of time series equals g*omega/(2*pi). g >=2.
% boundary_type 1 x 1. 1:  boundary is like surface of ocean where there are no pressure fluctuations.
%                          There is a direct path and a path reflecting from surface 
%                      2: boundary is hard. Reflected path bounces from hard flat rocks on bottom of ocean.
% fign          1 x 1. Figure number for plots. If empty, no plots shown
%
% OUTPUTS
%
% h             Structured array with fields
%
%   -taxis      1 x n.  t(i) is i'th time of data at receiver (s)
%   -s_emit     1 x n.  s_emit(i) is emitted signal at source at time  t(i).  This is the signal from the source
%               in absence of any other signals.
%   -s          1 x n.  s(i) is pressure fluctuation of signal at receiver. s(i) is at  time taxis(i).
%   -c3d        1 x 1. 3D effective speed (m/s) measured from time of emitted  peak
%               to time of  biggest received peak. Equals l_d/(t_r-t_emit)
%   -t_emit     1 x 1. Time (s) when emitted signal reaches peak
%   -t_r        1 x 1. Time (s) when received signal reaches maximum
%   -l_direct   1 x 1. Pythagorean distance (m)  between source and receiver.
% ierr          0: No error detected. 1: otherwise
%
function [h,ierr]=exact_time_domain_solution_reflected_direct_paths(c,src_xz,rec_xz,tau,omega,g,boundary_type,fign)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preliminary Steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
progname='exact_time_domain_solution_reflected_direct_paths.m';
% fprintf(1,'%s: See jlsy.mat to  debug. Pausing\n',progname);save  jlsy;pause

dotest=0;
if(dotest==1)
    fprintf(1,'%s: Test case. Press Enter to continue\n',progname);pause
    c=1500;
    
    scenario=8;
    g=10000;
    switch scenario
        case 1
            % Surface reflecting path
            % This yields c3d=1341.697620 m/s
            boundary_type=1;
            src_xz=[0 -40];
            rec_xz=[1 -5];
            tau=0.025;
            omega=100*2*pi;
        case 2
            % Surface reflecting path
            % This yields c3d=1229.m/s
            boundary_type=1;
            src_xz=[0 -20];
            rec_xz=[1 -5];
            tau=0.025;
            omega=100*2*pi;
        case 3
            % Surface reflecting path
            % This yields c3d=1093 m/s
            boundary_type=1;
            src_xz=[0 -10];
            rec_xz=[1 -5];
            tau=0.025;
            omega=100*2*pi;
        case 4
            % Surface reflecting path
            % This yields c3d=1226 m/s
            boundary_type=1;
            src_xz=[0 -2];
            rec_xz=[1 -5];
            tau=0.025;
            omega=100*2*pi;
            
        case 5
            % Surface reflecting path
            % This yields c3d=896 m/s
            boundary_type=1;
            src_xz=[0 -4.9];
            rec_xz=[1 -5];
            tau=0.025;
            omega=100*2*pi;
            
        case 6
            % Surface reflecting path
            % This yields c3d=1567 m/s
            boundary_type=1;
            src_xz=[0 -1.9];
            rec_xz=[1 -2];
            tau=0.025;
            omega=100*2*pi;
            
        case 7
            % Surface reflecting path
            % This yields c3d=980 m/s
            boundary_type=1;
            src_xz=[0 -4.9];
            rec_xz=[0 -5];
            tau=0.025;
            omega=100*2*pi;
            
        case 8
            % Surface reflecting path
            % This yields c3d=1807 m/s
            boundary_type=1;
            src_xz=[0 -0.1];
            rec_xz=[0 -5];
            tau=0.025;
            omega=100*2*pi;
        case 9
            % Surface reflecting path
            % This yields c3d=8596 m/s
            boundary_type=1;
            src_xz=[0 -0.1];
            rec_xz=[0 -5];
            tau=0.025;
            omega=20*2*pi;
        case 10
            % Surface reflecting path
            % This yields c3d=3792 m/s
            boundary_type=1;
            src_xz=[0 -4.9];
            rec_xz=[1 -5];
            tau=0.025;
            omega=20*2*pi;
            
            
        otherwise
            fprintf(1,'%s: software not written for scenario=%d\n',progname,scenario);
            pause
    end % switch scenario
    
    fign=1;
end

% Initialize outputs
h.taxis=[];
h.s_emit=[];
h.s=[];
h.c3d=[];
h.t_emit=[];
h.t_r=[];
h.l_direct=[];
ierr=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derive time bounds for emitted signal: t_bnds_emit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% multiplicative factor for tau
f=4;
t_bnds_emit=f*[-tau tau];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derive propagation time for reflected path to reach receiver: t_refl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get length of reflected path: l_refl
l_refl=sqrt((src_xz(1)-rec_xz(1))^2 + (src_xz(2) + rec_xz(2))^2 );
% Derive propagation time of reflected path: t_r
t_refl=l_refl/c;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get time  bounds for signal at receiver: t_bnds_rec
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_bnds_rec=[t_bnds_emit(1) t_refl+f*tau];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derive time axis for all signals: taxis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derive sample frequency: fs
% g multiplies highest frequency. Must be >=2
% g=10000;
fs=g*omega/(2*pi);
% Derive sample interval
dt=1/fs;

taxis=t_bnds_rec(1):dt:t_bnds_rec(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate emitted signal for all values in taxis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h.s_emit=exp(-(taxis/tau).^2).*cos(omega*taxis);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate received time series  from direct path: s_direct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get distance of straight path from source to receiver: l_direct
l_direct=sqrt((src_xz(1) - rec_xz(1))^2 + (src_xz(2) - rec_xz(2))^2);
% Get delay of direct path: t_delay_direct
t_delay_direct=l_direct/c;

% Get amplitude due to spherical  spreading of energy: a_direct
a_direct=1/l_direct;
% Derive signal at receiver from direct  path
s_direct=a_direct * exp(-((taxis-t_delay_direct)/tau).^2).*cos(omega*(taxis-t_delay_direct));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate received time series  from reflected path: s_reflect
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get amplitude of reflected path
a_refl=1/l_refl;
if(boundary_type==1)
    % Path reflections from pressure-release surface
    sign_use=-1;
else
    sign_use=1;
end

s_reflect = sign_use * a_refl*exp(-((taxis-t_refl)/tau).^2).*cos(omega*(taxis-t_refl));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate received time series 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s=s_direct + s_reflect;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate envelopes of emitted and received signals: ah_s_emit & ah_s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get envelope of emitted signal
ah_s_emit=abs(hilbert(h.s_emit));
%  Get envelope of received signal
ah_s=abs(hilbert(s));
% Get index of time for peak of emitted signal
[~,ind_s_emit_max]=max(ah_s_emit);
% Get index of time for peak of received signal
[~,ind_s_max]=max(ah_s);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute c3d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get time emitted peak leaves source
t_emit=taxis(ind_s_emit_max);
% Get time  received signal reaches peak
t_r=taxis(ind_s_max);
% Compute time for signal peak at receiver
t_src_rec_pk=t_r-t_emit;

%--------------------------------------------------------------------------------------------------------------
% Check direct path theoretical time, t_delay_direct, is very close to peak of envelope of just direct
% signal if there is no reflected path
%--------------------------------------------------------------------------------------------------------------
% Get envelope of signal at receiver if ONLY direct path is received
ah_s_direct=abs(hilbert(s_direct));
% Get index of its peak
[~,ind_ah_s_direct_peak]=max(ah_s_direct);
% Get time of this peak
t_ah_s_direct_peak=taxis(ind_ah_s_direct_peak);
% Derive increment of taxis
% dt_taxis=diff(taxis(1:2));
% Get absolute value of difference between theoretical and measured times of direct path
dt_direct=abs(t_ah_s_direct_peak - t_delay_direct);
% Get fractional error of dt_direct
fract_t_direct=dt_direct/t_delay_direct;
if(fract_t_direct > 0.001)
    fprintf(1,'%s: Measured time of direct path to receiver without interference is\n',progname);
    fprintf(1,' too far away from its theoretical value.  Probably need to decrease\n');
    fprintf(1,' sample interval of time axis by increasing g\n');
    fprintf(1,' diff between theoretical and measured peak = %e s.\n',abs(t_ah_s_direct_peak - t_delay_direct));
    fprintf(1,' Ratio of this difference to theoretical time direct path = %f\n',fract_t_direct);
    ierr=1;
    return
end
% Compute 3D effective speed of sound
c3d=l_direct/t_src_rec_pk;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h.taxis=taxis;
h.s=s;
h.c3d=c3d;
h.t_emit=taxis(ind_s_emit_max);
h.t_r=taxis(ind_s_max);
h.l_direct=l_direct;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot if  wanted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(~isempty(fign))
    figure(fign);
    clf;
    %-------------------------------------------------------------------------------------------------------------------
    % Plot emitted signal
    %-------------------------------------------------------------------------------------------------------------------
    ax1=subplot(4,1,1);
    plot(taxis,h.s_emit,'k--');
    hold on
    % Plot envelope and its max value
    plot(taxis,ah_s_emit,'k-');
    plot(taxis(ind_s_emit_max),ah_s_emit(ind_s_emit_max),'r.');

    xlabel('(s)');
    title('Emitted (dash) Envelope (solid)');
   
    %-------------------------------------------------------------------------------------------------------------------
    % Plot received signal
    %-------------------------------------------------------------------------------------------------------------------
    ax4=subplot(4,1,2);
    plot(taxis,s,'k--');
    xlabel('(s)');
    hold on;
    % plot its envelope and peak time
    plot(taxis,ah_s,'k-');
    plot(taxis(ind_s_max),ah_s(ind_s_max),'r.');
    title('Received Signal');
    
    %-------------------------------------------------------------------------------------------------------------------
    % Plot direct signal
    %-------------------------------------------------------------------------------------------------------------------
    ax2=subplot(4,1,3);
    plot(taxis,s_direct,'k--');
    hold on
   
    plot(taxis,ah_s_direct,'k-');
    % Get peak of ah_s_direct
    [pk_ah_s_direct,kk_direct]=max(ah_s_direct);
    plot(taxis(kk_direct),pk_ah_s_direct,'r.');
    xlabel('(s)');
    title('Direct Path at Receiver');
    
    %-------------------------------------------------------------------------------------------------------------------
    % Plot reflected signal
    %-------------------------------------------------------------------------------------------------------------------
    ax3=subplot(4,1,4);
    plot(taxis,s_reflect,'k--');
    hold on
    xlabel('(s)');
    title('Reflected path at Receiver');
        % Get envelope
    ah_s_reflect=abs(hilbert(s_reflect));
    plot(taxis,ah_s_reflect,'k-');
    % Get peak of ah_s_reflect
    [pk_ah_s_reflect,kk_reflect]=max(ah_s_reflect);
    plot(taxis(kk_reflect),pk_ah_s_reflect,'r.');
    
    

    linkaxes([ax1, ax2, ax3, ax4],'x');
    fprintf(1,'c3d=%f m/s\n',c3d);
    
    fprintf(1,'Peak time of envelope emitted signal= %f s\n',t_emit);
    fprintf(1,'Peak time of envelope received signal= %f s\n',t_r);
    fprintf(1,'Time for direct path to reach receiver = %f s\n',t_delay_direct);
    fprintf(1,'Time of peak of direct path envelope=%f s\n',taxis(kk_direct));
    fprintf(1,'Time for surf refl path to reach receiver= %f s\n',t_refl);
    fprintf(1,'Time of peak of reflected path envelope=%f s\n',taxis(kk_reflect));
    fprintf(1,' Amplitudes of direct and reflected paths = %f %f\n',a_direct,a_refl);
    
end

