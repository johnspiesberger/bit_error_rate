% Simulates two waveforms, each representing a specific symbol. Follow ideas in Stenner, Gauthier, and Neifeld
% letters to nature, vol 425, p. 695-8, 2003 = paper1.
%
% INPUTS
%
% f               1 x 1. Center frequency (Hz).  Can be zero if no carrier frequency desired.
% tau             1 x 1. Time scale of initial Gaussian envelope before and up to its peak (s).
%                 Gaussian envelope is exp(-(t/tau)^2) so peak of envelope occurs at t=0. See Fig 2a, paper1
%                 before point of non-a
% tau_after       1 x 3. Same as tau except pertains to Gaussian envelopes at t>=0. See. Fig. 2a, paper1.  t=0
%                 corresponds to the point of non-analyticity.
%                 tau_after(1) pertains to waveform (symbol 0) that goes down at t=0. tau_after(2) to waveform 
%                 (symbol 1) that goes up at t=0 to its peak.  tau_after(3) pertains to symbol 1 after it reaches its
%                 peak.
% amp_big         1 x 1.  Maximum amplitude of large amplitude waveform after point of non-analyticity.
%                 Is > 1.  amp_big corresponds to about 1.5 in Fig. 2a of paper1.
%                 carrier is maximum where the Gaussian is a maximum.  
%                 The waveform before t=0 plus this Gaussian is symbol 1.
%                 Symbol "0"  corresponds to the Gaussian envelope that goes down from time 0.
% f_samp          1 x 1. Sample frequency (Hz) of output waveforms.  Must be > 2*f
% n_tau_to_zero   1 x 1. Gaussian envelope is set to zero at t=-n_tau_to_zero * tau for initial part of waveform
%                 and is zero at delay n_tau_to_zero * tau_after after peak is reached for each waveform.
%                 Typical value might be 3, thus exp(-(3*tau/tau)^2) = exp(-9) = 1.2d-4.
% fign_ideal_symb 1 x 1. Matlab figure number to plot ideal emitted symbols. If empty, plot not made
%
% OUTPUTS
%
% w               2 x q.  w(i,:) is time series of ideal emitted signal i.  This does  not include bandpass filter
%                 to mimic real finite bandwidth effects of any transmitter. Sample frequency is f_samp.
% e               2 x q.  Same as w except it the Gaussian envelope only of the ideal emitted pulse.
% t               1 x q.  w(i,j) occurs at time t(j).  t(j)=0 corresponds to point of non-analyticity: where the
%                 two ideal emitted pulses are switched from symbol "0" or symbol "1".  
% ierr            0: no errors detected
%                 1: otherwise
%
function [w,e,t,ierr]=simulate_two_waveforms(f,tau,tau_after,amp_big,f_samp,n_tau_to_zero,fign_ideal_symb)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preliminary Steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
progname='simulate_two_waveforms.m';
% fprintf(1,'%s: See jls7u.mat to debug. pausing\n',progname);pause

dotest=0;
if(dotest==1)
    fprintf(1,'%s: Test case. Press ENTER to continue\n',progname);pause
    f=0;
    tau=0.025;
    tau_after=[0.01 0.03 .05];
    amp_big=1.5;
    f_samp=1.d4; 
    n_tau_to_zero=3;
end

% Intialize outputs
w=[];
e=[];
t=[];
ierr=0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct waveform up to time zero: w0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get approximate first time to simulate: t_first
t_first_sugg=-n_tau_to_zero * tau;
% Get number samples between t_first and t=0
n=ceil(abs(t_first_sugg)*f_samp);
% Get first sample time
t_first=-(n-1)/f_samp;

% Get time axis between t_first and t=0
t0=linspace(t_first,0,n);
% Get waveform
w0=exp(-(t0/tau).^2) .* cos(2*pi*f*t0);
% Get its envelope
e0=exp(-(t0/tau).^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct waveform after t=0 for symbol 1 which goes up in amplitude: w_symb_1 & t_after
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------------------------------------------------------------------------------------------------------
% Find time at which symbol 1 reaches max amplitude of amp_big: t_max
%-----------------------------------------------------------------------------------------------------------------------
% Solve 
% amp_big*exp(-(t_max/tau_after)^2) = 1
% 1/amp_big=exp(-(t_max/tau_after)^2)
% ln(1/amp_big) = -(t_max/tau_after)^2
% sqrt(-ln(1/amp_big)) = t_max/tau_after
% t_max=tau_after * sqrt(-ln(1/amp_big))            (1)
%
t_max=tau_after(2) * sqrt(-log(1/amp_big));

%-----------------------------------------------------------------------------------------------------------------------
% Find time interval where symbol 1 goes from its peak down to exp(-n_tau_to_zero^2): dt_last
%-----------------------------------------------------------------------------------------------------------------------
% Solve this for dt_last   
% amp_big*exp(-(dt_last/tau_after)^2) = exp(-n_tau_to_zero^2)
%  exp(-(dt_last/tau_after)^2) = exp(-n_tau_to_zero^2) / amp_big
% -(dt_last/tau_after)^2 = log(exp(-n_tau_to_zero^2) / amp_big)
%  dt_last = sqrt( -log(exp(-n_tau_to_zero^2) / amp_big)) * tau_after    (2)

dt_last = sqrt( -log(exp(-n_tau_to_zero^2) / amp_big)) * tau_after(3);

%-----------------------------------------------------------------------------------------------------------------------
% Derive last time for symbol 1: t_to_end
% Derive time axis for t>0: t_after
%-----------------------------------------------------------------------------------------------------------------------
% Get interval between t=0 and end of time series: t_to_end
t_to_end=t_max + dt_last;
% Get approximate number samples from t=0 to t_to_end
nn=ceil(t_to_end*f_samp);
% Get time axis after t=0 to t_to_end
t_after=1/f_samp + (0:nn-1)/f_samp;
%-----------------------------------------------------------------------------------------------------------------------
% Construct time series for symbol 1 after t=0: w_symb_1
%-----------------------------------------------------------------------------------------------------------------------
% Find indices of t_after from t=0 to peak of symbol:  t_after(p1)
[p1]=find(t_after<=t_max);
% Find indices of t_after  corresponding to > t_max: t_after(p2)
[p2]=find(t_after>t_max);

w_symb_1 = [...
    amp_big*exp(-((t_after(p1)-t_max)/tau_after(2)).^2) .* cos(2*pi*f*t_after(p1)) ...
    amp_big*exp(-((t_after(p2)-t_max)/tau_after(3)).^2) .* cos(2*pi*f*t_after(p2))];

% Get Gaussian envelope of symbol
e_symb_1= [...
    amp_big*exp(-((t_after(p1)-t_max)/tau_after(2)).^2) ...
    amp_big*exp(-((t_after(p2)-t_max)/tau_after(3)).^2)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct waveform after t=0 for symbol 0 which goes down in amplitude: w_symb_0 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w_symb_0 = exp(-(t_after/tau_after(1)).^2) .* cos(2*pi*f*t_after);

% Get Gaussian envelope of symbol
e_symb_0= exp(-(t_after/tau_after(1)).^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct symbols 0 & 1 for entire time: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get time axis for all symbols
t=[t0 t_after];
% Get # time samples
q=length(t);
% Allocate memory for output symbols
w=nan(2,q);
% Fill it
w(1,:)=[w0 w_symb_0];
w(2,:)=[w0 w_symb_1];

%-----------------------------------------------------------------------------------------------------------------------
% Construct envelopes for entire time
%-----------------------------------------------------------------------------------------------------------------------
% Allocate memory for output
e=nan(2,q);
% Fill it
e(1,:)=[e0 e_symb_0];
e(2,:)=[e0 e_symb_1];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot both symbols
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(~isempty(fign_ideal_symb))
    figure(fign_ideal_symb);
    clf;
    subplot(2,1,1);
    plot(t,w(1,:),'g-',t,w(2,:),'r-');

    xlabel('(s)');
    ylabel('Output');
    title('Ideal emitted symbols. Envelopes (dash)');
        
    subplot(2,1,2);
    plot(t,e(1,:),'g-',t,e(2,:),'r-');
    xlabel('(s)');
    ylabel('Output');
    title('Ideal emitted envelope of symbols. Envelopes (dash)');
    
    
end

