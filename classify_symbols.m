% Classifier for received symbols 0 and 1 in one of two conditions: 1) no boundary reflection. Direct path only.
% 2) boundary reflection.  direct and reflected paths interfere. This emperically tests bit_error_rate.m
% Ideas follow the procedure outlined in Fig. 3 of
% Stenner, Gauthier, and Neifeld letters to nature, vol 425, p. 695-8, 2003 (paper 1).
%
% INPUTS
%
% symbol_0_no_noise    1 x nt. Time series at receiver due to symbol 0 without noise. Includes all effects, including
%                      bandpass due to optical switch and the receiver's bandwidth (if receiver's bandwidth is modeled).
% symbol_1_no_noise    1 x nt. Same as symbol_0_no_noise except is symbol 1.
% dt                   1 x 1. Time interval (s) between tabulated values in symbol_0 & symbol_1
% t_s                  1 x 1. (s). Defined in Fig. 3 of paper 1 above.
% ind_direct_na_arrive 1 x 1. Suppose point of non-analyticity arrives at receiver along direct path at indices
%                      ind0 & ind1.  Then ind_direct_na_arrive=min(ind0,ind1).
%                      
% ind_tau_alpha        1 x 1. Corresponds to index of symbol_0_no_noise & symbol_1_no_noise before point where 
%                      these symbols separate in time. Output bit error rate is tabulated starting at this index.
% ind_na_arrive_no_filters_no_interference
%                      1 x 1. Index of symbol_0_no_noise and symbol_1_no_noise corresponding to time when
%                      point of non-analyticity would arrive at receiver if said point traveled to receiver
%                      at phase speed of signal without any effects due to temporal interference of finite bandwidth
%                      filtering.   I.e. it corresponds to time when point of non-analyticity of electromagnetic wave 
%                      would arrive if it propagated at the speed of light in a vacuum without any interference
%                      effects.   If symbol 0 or 1 can be detected at time before this index, there is theoretical 
%                      evidence that special relaitivity might be violated and causility is invalidated.
% snr                  1 x 1. Signal-to-noise ratio (dB) at receiver for symbol 1 at its peak value.
% n_symbol_transmit    1 x 1. Defined by simulate_c_information_measurement.m. In this program, equals the
%                      number of times each noisy symbol is simulated for purpose of classifying received signal.
% n_samps_before_c
%                      1 x 1. Calculations of the bit error rate are tabulated starting at
%                      n_samps_before_c samples before the time the point of non-analyticity arrives at receiver
%                      in the absence of temporal interference, finite bandwidth filtering, & and noise.
%                      Is >=0.  Zero corresponds to time when light arrrives at receiver in vacuum when there is
%                      no interference (and not infinite bandwidth and no noise).  Want to make  n_samps_before_c
%                      large enough so you can see if detection of symbols can be discerened before light would get there
%                      in a vacuum. Sample interval is 1/f_samp (above).
% t_classify           1 x 1 Time (s) when decision is made for classification of symbol. This time corresponds to
%                      this many seconds past t_axis(ind_tau_alpha).
%
% t_axis               1 x nt. t_axis(i) is time (s) corresponding to symbol_0_no_noise(i). 
% fign                 1 x 1.  Matlab figure for plotting ber_info. If emtpy, no plot made
%
% OUTPUTS
%
% ber_emperical        Structured array with  fields
%   -ber_symb_0        1 x 1. Bit error rate when noisy symbol  0 is present 
%   -ber_symb_1        1 x 1. Bit error rate when noisy symbol  1 is present
%   -ber               1 x 1. Bit error rate
% ierr                 0: no error detected
%                      1; otherwise
%
function [ber_emperical,ierr]=classify_symbols(symbol_0_no_noise,symbol_1_no_noise,dt,t_s,ind_direct_na_arrive,...
    ind_tau_alpha,ind_na_arrive_no_filters_no_interference,snr,n_symbol_transmit,n_samps_before_c,t_classify,t_axis,...
    fign)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preliminary Steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
progname='classify_symbols.m';
% fprintf(1,'%s: See jls7ub.mat to debug. pausing\n',progname);save jls7ub;pause


fprintf(1,'%s: This program gave  same answer as bit_error_rate.m changeset 14 when noise\n',progname);
fprintf(1,' for white noise. This routine could be further modified to filter the noise as well.\n');
fprintf(1,' Press ENTER  to continue\n');


% Initialize outputs
ber_emperical.ber_symb_0=[];
ber_emperical.ber_symb_1=[];
ber_emperical.ber=[];
ierr=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derive index of t_axis where  classification of symbol is made: ind_classify
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code from ber_error_rate.m computes value of ber_info.tau and ber starting with t_axis(ind_first).
% ind_first=ind_na_arrive_no_filters_no_interference - n_samps_before_c;

% In bit_error_rate.m, n_samps_before_c + 1 is the time axis index corresponding to 0 ber_info.tau=0
% ber_error_rate.m states  ber_info.tau=0 occurs at t_axis(ind_tau_alpha) 
[~,ind_classify]=min(abs(t_axis - t_classify - t_axis(ind_tau_alpha)));

% Code below identical to bit_error_rate.m except as noted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BER for direct path mimics scenario in vacuum
% "       symbol 0 (r0) mimics scenario of fast light propagation

% Get length of symbol_0_no_noise: nt
nt=length(symbol_0_no_noise);
%-----------------------------------------------------------------------------------------------------------------------
% Derive index of symbol_0_no_noise corresponding to t_s in paper 1: ind_t_s
%-----------------------------------------------------------------------------------------------------------------------
ind_t_s=round(t_s/dt) + 1;

%-----------------------------------------------------------------------------------------------------------------------
% Get standard deviation of noise to add to symbols without noise: std_n
%-----------------------------------------------------------------------------------------------------------------------
% Get max value of symbol 1 at receiver due to direct and reflected paths: a_max
a_max=max(symbol_1_no_noise);
% snr = 20log10(a_max/std_n)
% a_max=10^(snr/20) * std_n
std_n=a_max/10^(snr/20);

%-----------------------------------------------------------------------------------------------------------------------
% Compute N_j(tau); j=0,1.  Defined in paper 1
%-----------------------------------------------------------------------------------------------------------------------
% N_0(tau) integrates symbol_0_no_noise_0^2 because it is symbol 0
N_0=cumtrapz(symbol_0_no_noise(ind_t_s:nt).^2) * dt;
% N_1(tau) integrates symbol_1_no_noise^2 because it is symbol 1
N_1=cumtrapz(symbol_1_no_noise(ind_t_s:nt).^2) * dt;

%-----------------------------------------------------------------------------------------------------------------------
% Get N_j(tau_alpha); j=0,1.  Defined in paper 1
%-----------------------------------------------------------------------------------------------------------------------
N_0_tau_alpha=N_0(ind_tau_alpha);
N_1_tau_alpha=N_1(ind_tau_alpha);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute alpha_j(tau_alpha); j=-0,1. paper1:  alpha_0_tau_alpha & alpha_1_tau_alpha
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get number indices in [ind_t_s , ind_tau_alpha]: nn
nn=ind_tau_alpha - ind_t_s+1;

% Initialize
alpha_0_tau_alpha=0;
alpha_1_tau_alpha=0;
% Derive alpha_0_tau_alpha & alpha_1_tau_alpha by using mean values over n_symbol_transmit realizations of noise
for i=1:n_symbol_transmit
    % Generate noisy data for symbol_0_no_noise & symbol_1_no_noise: chi)0 & chi_1
    chi_0=symbol_0_no_noise(ind_t_s:ind_tau_alpha) + std_n*randn(1,nn);
    chi_1=symbol_1_no_noise(ind_t_s:ind_tau_alpha) + std_n*randn(1,nn);
    % Accumulate integral for i'th case of noise
    alpha_0_tau_alpha = alpha_0_tau_alpha + trapz(chi_0.^2) * dt;
    alpha_1_tau_alpha = alpha_1_tau_alpha + trapz(chi_1.^2) * dt;
end % for i=1:n_symbol_transmit

alpha_0_tau_alpha=alpha_0_tau_alpha/n_symbol_transmit / N_0_tau_alpha;
alpha_1_tau_alpha=alpha_1_tau_alpha/n_symbol_transmit / N_1_tau_alpha;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine number values of tau for: n_tau
n_tau=nt-ind_tau_alpha + 1;


% Get first index of symbols to compute BER for
ind_first=ind_na_arrive_no_filters_no_interference - n_samps_before_c;

% Get # times to compute D for: nD
% nD=nt-ind_tau_alpha+1;
nD=nt-ind_first + 1;


% Allocate memory to store D(tau) = l_0(tau) - l_1(tau) for tau corresponding to ind_tau_alpha:nt 
% D0 uses matched filter with replica symbol_0_no_noise
% D1 "                                symnbol_1_no_noise
D0=zeros(n_symbol_transmit,nD,'single');
D1=zeros(n_symbol_transmit,nD,'single');

% Allocate memory to store D0 and D1 another way for parfor
% g.sim=[];
% D0=repmat(g,1,nD);
% D1=repmat(g,1,nD);
% clear g


for i=1:n_symbol_transmit
    if(rem(i,10)==0)
        fprintf(1,'Starting simulation %d of %d\n',i,n_symbol_transmit);
    end
    % Generate noisy data for symbol_0_no_noise & symbol_1_no_noise: chi)0 & chi_1
    chi_0=symbol_0_no_noise + std_n*randn(1,nt);
    chi_1=symbol_1_no_noise + std_n*randn(1,nt);
    % chi_0=symbol_0_no_noise;
    % chi_1=symbol_1_no_noise;

    
    % Loop over ind_now
    % for ind_now=ind_tau_alpha:nt
    % Set ind_now to ind_classify for classify_symbols.m and comment out for ind_nnow=ind_first:nt
    %for ind_now=ind_first:nt
        ind_now=ind_classify;
        %---------------------------------------------------------------------------------------------------------------
        % Get D  when only symbol_0_no_noise + noise    are present
        %---------------------------------------------------------------------------------------------------------------
        % Replicas for matched filter are symbol_0_no_noise and symbol_1_no_noise
        % Derive realization of l_0 for ind_t_s:ind_now:  this is a matched filter
        l_0=trapz(chi_0(ind_t_s:ind_now).*symbol_0_no_noise(ind_t_s:ind_now))/(alpha_0_tau_alpha*N_0(ind_now));
        % Derive realization of l_1 for ind_t_s:ind_now:  this is a matched filter
        l_1=trapz(chi_0(ind_t_s:ind_now).*symbol_1_no_noise(ind_t_s:ind_now))/(alpha_1_tau_alpha*N_1(ind_now));
        
        % Get index of D to store ind_now: ind_D
        % ind_D=ind_now-ind_tau_alpha+1;
        ind_D=ind_now-ind_first+1;
        D0(i,ind_D)=single(l_0-l_1);
        % D0(i).sim=single(l_0-l_1);
        %---------------------------------------------------------------------------------------------------------------
        % Get D when only symbol_1_no_noise + noise are present
        %---------------------------------------------------------------------------------------------------------------
        % Replicas for matched filter are symbol_0_no_noise and symbol_1_no_noise
        % Derive realization of l_0 for ind_t_s:ind_now:  this is a matched filter
        l_0=trapz(chi_1(ind_t_s:ind_now).*symbol_0_no_noise(ind_t_s:ind_now))/(alpha_0_tau_alpha*N_0(ind_now));
        % Derive realization of l_1 for ind_t_s:ind_now:  this is a matched filter
        l_1=trapz(chi_1(ind_t_s:ind_now).*symbol_1_no_noise(ind_t_s:ind_now))/(alpha_1_tau_alpha*N_1(ind_now));
        
        % Get index of D to store ind_now: ind_D
        % ind_D=ind_now-ind_tau_alpha+1;
        
        D1(i,ind_D)=single(l_0-l_1);
      
    % end % for ind_now=ind_tau_alpha:nt
    
end % for i=1:n_symbol_transmit

%-----------------------------------------------------------------------------------------------------------------------
% Code below is not in bit_error_rate.m
%-----------------------------------------------------------------------------------------------------------------------
% use_normalizaition: TRUE Use normalization of matched filters as  in Stenner et al
%                     False: no normalization of matched filters
use_normalization=true(1,1);
% Set state of random number generator to see differences with normalization and not normalization.
rng(1);
fprintf(1,'%s: Set rng(1). Press ENTER to continue.\n',progname);pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute output of matched filter at t_classify when input symbol has no noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind_now=ind_classify;
if(use_normalization)
    %-----------------------------------------------------------------------------------------------------------------------
    % Get outputs of both matched filters when symbol 0 inputted
    %-----------------------------------------------------------------------------------------------------------------------
    % matched filter for symbol 0 no noise
    mf_perfect_symb_0_in_symbol_0_replica=...
        trapz(symbol_0_no_noise(ind_t_s:ind_now).*symbol_0_no_noise(ind_t_s:ind_now))/(alpha_0_tau_alpha*N_0(ind_now));
    % matched filter for symbol 1 no noise
    mf_perfect_symb_0_in_symbol_1_replica=...
        trapz(symbol_0_no_noise(ind_t_s:ind_now).*symbol_1_no_noise(ind_t_s:ind_now))/(alpha_1_tau_alpha*N_1(ind_now));
    %-----------------------------------------------------------------------------------------------------------------------
    % Get outputs of both matched filters when symbol 1 inputted
    %-----------------------------------------------------------------------------------------------------------------------
    mf_perfect_symb_1_in_symbol_0_replica=...
        trapz(symbol_1_no_noise(ind_t_s:ind_now).*symbol_0_no_noise(ind_t_s:ind_now))/(alpha_0_tau_alpha*N_0(ind_now));
    % Derive realization of l_1 for ind_t_s:ind_now:  this is a matched filter
    mf_perfect_symb_1_in_symbol_1_replica=...
        trapz(symbol_1_no_noise(ind_t_s:ind_now).*symbol_1_no_noise(ind_t_s:ind_now))/(alpha_1_tau_alpha*N_1(ind_now));
else
    % No normalization
    %-----------------------------------------------------------------------------------------------------------------------
    % Get outputs of both matched filters when symbol 0 inputted
    %-----------------------------------------------------------------------------------------------------------------------
    % matched filter for symbol 0 no noiise
    mf_perfect_symb_0_in_symbol_0_replica=...
        trapz(symbol_0_no_noise(ind_t_s:ind_now).*symbol_0_no_noise(ind_t_s:ind_now));
    % matched filter for symbol 1 no noise
    mf_perfect_symb_0_in_symbol_1_replica=...
        trapz(symbol_0_no_noise(ind_t_s:ind_now).*symbol_1_no_noise(ind_t_s:ind_now));
    %-----------------------------------------------------------------------------------------------------------------------
    % Get outputs of both matched filters when symbol 1 inputted
    %-----------------------------------------------------------------------------------------------------------------------
    mf_perfect_symb_1_in_symbol_0_replica=...
        trapz(symbol_1_no_noise(ind_t_s:ind_now).*symbol_0_no_noise(ind_t_s:ind_now));
    % Derive realization of l_1 for ind_t_s:ind_now:  this is a matched filter
    mf_perfect_symb_1_in_symbol_1_replica=...
        trapz(symbol_1_no_noise(ind_t_s:ind_now).*symbol_1_no_noise(ind_t_s:ind_now));

end % if(use_normalization)

% Print values
fprintf(1,'\n');
fprintf(1,'Noiseless symbol inputs to matched_filter\n');
fprintf(1,'Matched filt symb 0 with symb 0 replica=%e\n',mf_perfect_symb_0_in_symbol_0_replica);
fprintf(1,'Matched filt symb 0 with symb 1 replica=%e\n',mf_perfect_symb_0_in_symbol_1_replica);
fprintf(1,'\n');
fprintf(1,'Matched filt symb 1 with symb 0 replica=%e\n',mf_perfect_symb_1_in_symbol_0_replica);
fprintf(1,'Matched filt symb 1 with symb 1 replica=%e\n',mf_perfect_symb_1_in_symbol_1_replica);
% It may be the perfect matched filter is larger for
% symb 0 with  replica 1 than symb 1 with  replica 1 OR
% "    1               0 than      1 with  replica 1
% because the whole symbol is not yet received.
D0_perfect=mf_perfect_symb_0_in_symbol_0_replica - mf_perfect_symb_0_in_symbol_1_replica;

D1_perfect=mf_perfect_symb_1_in_symbol_0_replica - mf_perfect_symb_1_in_symbol_1_replica;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute output of matched filter at t_classify with noisy symbols
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allocate memory to store outputs of four matched filters
symb_0_in_symbol_0_replica=nan(1,n_symbol_transmit);
symb_0_in_symbol_1_replica=nan(1,n_symbol_transmit);
symb_1_in_symbol_0_replica=nan(1,n_symbol_transmit);
symb_1_in_symbol_1_replica=nan(1,n_symbol_transmit);

% Allocate memory to store output of matched filter for symbol 0 due only to noise
noise_output_symbol_0_replica=nan(1,n_symbol_transmit);

for i=1:n_symbol_transmit

    % Generate noise for symbol 0
    noise_symb_0=std_n*randn(1,nt);
    % Generate noise for symbol 1
    noise_symb_1=std_n*randn(1,nt);

    % Generate noisy data for symbol_0_no_noise & symbol_1_no_noise: chi_0 & chi_1
    chi_0=symbol_0_no_noise + noise_symb_0;
    chi_1=symbol_1_no_noise + noise_symb_1;

    if(use_normalization)
        % Apply matched filters to noisy symbol 0
        symb_0_in_symbol_0_replica(i)  =trapz(chi_0(ind_t_s:ind_now).*symbol_0_no_noise(ind_t_s:ind_now))/(alpha_0_tau_alpha*N_0(ind_now));
        symb_0_in_symbol_1_replica(i)=trapz(chi_0(ind_t_s:ind_now).*symbol_1_no_noise(ind_t_s:ind_now))/(alpha_1_tau_alpha*N_1(ind_now));
        % Apply matched filters to noisy symbol 1
        symb_1_in_symbol_0_replica(i)  =trapz(chi_1(ind_t_s:ind_now).*symbol_0_no_noise(ind_t_s:ind_now))/(alpha_0_tau_alpha*N_0(ind_now));
        symb_1_in_symbol_1_replica(i)=trapz(chi_1(ind_t_s:ind_now).*symbol_1_no_noise(ind_t_s:ind_now))/(alpha_1_tau_alpha*N_1(ind_now));

        % Compute noise only output of matched filter for symbol 0
        noise_output_symbol_0_replica(i) = trapz(noise_symb_0(ind_t_s:ind_now).*symbol_0_no_noise(ind_t_s:ind_now))/(alpha_0_tau_alpha*N_0(ind_now));
    else
        % No normalization

        % Apply matched filters to noisy symbol 0
        symb_0_in_symbol_0_replica(i)  =trapz(chi_0(ind_t_s:ind_now).*symbol_0_no_noise(ind_t_s:ind_now));
        symb_0_in_symbol_1_replica(i)=trapz(chi_0(ind_t_s:ind_now).*symbol_1_no_noise(ind_t_s:ind_now));
        % Apply matched filters to noisy symbol 1
        symb_1_in_symbol_0_replica(i)  =trapz(chi_1(ind_t_s:ind_now).*symbol_0_no_noise(ind_t_s:ind_now));
        symb_1_in_symbol_1_replica(i)=trapz(chi_1(ind_t_s:ind_now).*symbol_1_no_noise(ind_t_s:ind_now));

        % Compute noise only output of matched filter for symbol 0
        noise_output_symbol_0_replica(i) = trapz(noise_symb_0(ind_t_s:ind_now).*symbol_0_no_noise(ind_t_s:ind_now));

    end

end % for i=1:n_symbol_transmit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get BER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Use D as decision statistic.  For each noisy symbol, compute D, and classify the symbol as a 0 if
%                 closer to D0 than D1, and classify as a 1 if closer to D1 than  D0 where
%                 D0=mf_perfect_symb_0_in_symbol_0_replica - mf_perfect_symb_0_in_symbol_1_replica;
%                 D1=mf_perfect_symb_1_in_symbol_0_replica - mf_perfect_symb_1_in_symbol_1_replica;
%                 9Oct2025: Spiesberger thinks this is the correct method that mimics Stenner et al paper.
%
% Initialize number times D not closer to D0 than D1 when symbol 0 present
n_not_closer_D0_when_symb_0_present=0;
% Initialize number times D not closer to D1 than D0 when symbol 01 present
n_not_closer_D1_when_symb_1_present=0;

% Initialize symbol inputted to classifier.  0 means  symbol 0. 1  means symbol 1
symb_inputted_to_classifier=nan(1,n_symbol_transmit);
% Initialize memory to hold classification. 0 means classification is symbol  0
symbol_classification=nan(1,n_symbol_transmit);


% Fill them
for i=1:n_symbol_transmit
    % Compute D from noisy symbol 0
    symb_0_D=symb_0_in_symbol_0_replica(i) - symb_0_in_symbol_1_replica(i);
    % Compute D from noisy symbol 1
    symb_1_D=symb_1_in_symbol_0_replica(i) - symb_1_in_symbol_1_replica(i);

    % Determine if noisy symbol 0 correctly classified
    if(abs(symb_0_D-D0_perfect) > abs(symb_0_D-D1_perfect))
        % Noisy symbol 0 incorrectly classifed as symbol 1
        n_not_closer_D0_when_symb_0_present = n_not_closer_D0_when_symb_0_present + 1;
    else
        % Noisy symbol  0 correctly classified
    end

    % Determine if noisy symbol 1 correctly classified
    if(abs(symb_1_D-D0_perfect) <= abs(symb_1_D-D1_perfect))
        % Noisy symbol 1 incorrectly classifed as symbol 0
        n_not_closer_D1_when_symb_1_present = n_not_closer_D1_when_symb_1_present + 1;
    else
        % Noisy symbol  0 correctly classified
    end

    %------------------------------------------------------------------------------------------------------
    % Compute BER
    %------------------------------------------------------------------------------------------------------

    % Choose symbol 0 or 1 to use for input to D: symb_inputted_to_classifier(i)
    % Pretend I do not know what symbol is chosen.
    % There is a 50% chance of choosing symbol 0  assuming probability of sending 0 is 0.5.
    if(rand(1,1)<=0.5)
        symb_inputted_to_classifier(i)=0;
    else
        symb_inputted_to_classifier(i)=1;
    end

    % Get value of D for symbol inputted to classifier
    if(symb_inputted_to_classifier(i)==0)
        % D  is from noisy symbol 0
        D_use=symb_0_D;
    else
        % D is from noisy symbol 1
        D_use=symb_1_D;
    end

    % Classify symbol:
    if(abs(D_use-D0_perfect) <= (D_use-D1_perfect))
        % Classify as symbol 0
        symbol_classification(i)=0;
    else
        % Classify as symbol 1
        symbol_classification(i)=1;
    end


end % for i=1:n_symbol_transmit

% Get number times symbol incorrectly classified
n_bit_error = length(find(symbol_classification ~= symb_inputted_to_classifier));

% Transmitted n_symbol_transmit 0's and n_symbol_transmit 1's.
ber_emperical.ber_symb_0=n_not_closer_D0_when_symb_0_present/n_symbol_transmit;
ber_emperical.ber_symb_1=n_not_closer_D1_when_symb_1_present/n_symbol_transmit;

ber_emperical.ber = n_bit_error/n_symbol_transmit;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%9%%%%%%9%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-----------------------------------------------------------------------------------------------------------------------
% Plot histogram of noisy symbol 0 output from symbol  0 and symbol 1 matched filter
%-----------------------------------------------------------------------------------------------------------------------
V=zeros(4,4);
figure(50);
clf
% Plot output of symbol 0 mf with noisy symbol 0 input
subplot(2,1,1);
histogram(symb_0_in_symbol_0_replica);
title('symbol 0 mf with noisy symbol 0 input');
hold on;
V(1,:)=axis;
% Plot output of symbol 1 mf with noisy symbol 0 input
subplot(2,1,2);
histogram(symb_0_in_symbol_1_replica);
title('symbol 1 mf with noisy symbol 0 input');
hold on;
V(2,:)=axis;

% for i=1:2
%     subplot(2,1,i);
%     % Plot both  ideal matched filters for both symbols
%     plot(repmat(mf_perfect_symb_0_in_symbol_0_replica,1,2),V_use(3:4),'g-');
%         plot(repmat(mf_perfect_symb_0_in_symbol_1_replica,1,2),V_use(3:4),'r-');
%     xlabel('\color{green} replica 0 mf \color{red} replica 1 mf')
% 
% end

%-----------------------------------------------------------------------------------------------------------------------
% Plot histogram of noisy symbol 1 output from symbol  0 and symbol 1 matched filter
%-----------------------------------------------------------------------------------------------------------------------

figure(51);
clf
% Plot output of symbol 0 mf with noisy symbol 1 input
subplot(2,1,1);
histogram(symb_1_in_symbol_0_replica);
title('symbol 0 mf with noisy symbol 1 input');
hold on;
V(3,:)=axis;
% Plot output of symbol 1 mf with noisy symbol 1 input
subplot(2,1,2);
histogram(symb_1_in_symbol_1_replica);
title('symbol 1 mf with noisy symbol 1 input');
hold on;
V(4,:)=axis;

% for i=1:2
%     subplot(2,1,i);
%     % Plot both  ideal matched filters for both symbols
%     plot(repmat(mf_perfect_symb_1_in_symbol_0_replica,1,2),V_use(3:4),'g-');
%     plot(repmat(mf_perfect_symb_1_in_symbol_1_replica,1,2),V_use(3:4),'r-');
%     xlabel('\color{green} replica 0 mf \color{red} replica 1 mf')
% 
% end

% Get limits to use for both figures: V_use
V_use=[min(V(:,1)) max(V(:,2)) min(V(:,3)) max(V(:,4))];
% make both subplots have same axis limits for Figure 51
for i=1:2
    subplot(2,1,i);
    axis(V_use);
        % Plot both  ideal matched filters for both symbols
    plot(repmat(mf_perfect_symb_1_in_symbol_0_replica,1,2),V_use(3:4),'g-');
    plot(repmat(mf_perfect_symb_1_in_symbol_1_replica,1,2),V_use(3:4),'r-');
    xlabel('\color{green} replica 0 mf \color{red} replica 1 mf')

end

figure(50);
for i=1:2
    subplot(2,1,i);
    axis(V_use);
        plot(repmat(mf_perfect_symb_0_in_symbol_0_replica,1,2),V_use(3:4),'g-');
        plot(repmat(mf_perfect_symb_0_in_symbol_1_replica,1,2),V_use(3:4),'r-');
    xlabel('\color{green} replica 0 mf \color{red} replica 1 mf')


end
%-----------------------------------------------------------------------------------------------------------------------
% Plot D0 and D1
%-----------------------------------------------------------------------------------------------------------------------
% D0 and D1 show more discrimination for classification than  treating outputs of matched filter separately
D0=symb_0_in_symbol_0_replica - symb_0_in_symbol_1_replica;
D1=symb_1_in_symbol_0_replica - symb_1_in_symbol_1_replica;
% Allocate memory to store axes
V_all=zeros(2,4);
% Compute value of tau corresponding to j
% NO longer true:  tau=0 occurs at index ind_tau_alpha
% Now true: tau=0 corresponds to index ind_na_arrive_no_filters_no_interference
% tau_now=(j-1)*dt;
tau_now=( (ind_D - 1) - n_samps_before_c)*dt;

figure(100);
clf;
subplot(2,1,1);
hist0=histogram(D0,'normalization','probability');
title(['D only symbol 0 present. tau = ' num2str(tau_now) ' s. \color{red} perfect value']);
V_all(1,:)=axis;
hold on
subplot(2,1,2);
hist1=histogram(D1,'normalization','probability');
xlabel('\tau (s)');
title(['D only symbol 1 present. tau = ' num2str(tau_now) ' s. \color{red} perfect value']);
hold on
V_all(2,:)=axis;

% Make axis limits identical
V=[min(V_all(:,1)) max(V_all(:,2)) min(V_all(:,3)) max(V_all(:,4))];
axis(V);
subplot(2,1,1);
axis(V);

% Plot perfect results
subplot(2,1,1);
plot(repmat(D0_perfect,1,2),V_use(3:4),'r-');
subplot(2,1,2);
plot(repmat(D1_perfect,1,2),V_use(3:4),'r-');


%-----------------------------------------------------------------------------------------------------------------------
% Redraw D0 and D1 with identical bin edges
%-----------------------------------------------------------------------------------------------------------------------
%  Get Bin limits  to use for plots
BinLimits=[min([hist0.BinLimits hist1.BinLimits]) max([hist0.BinLimits hist1.BinLimits])];
% Get smallest bin interval
BinWidth=min([hist0.BinWidth hist1.BinWidth]);
% Get # bins to use for plot
NumBins=ceil(diff(BinLimits)/BinWidth);

figure(101);
clf;
% Plot D0
subplot(2,1,1);
hist00=histogram(D0,'normalization','probability','BinLimits',BinLimits,'BinWidth',BinWidth,'NumBins',NumBins);
title(['D only symbol 0 present. tau = ' num2str(tau_now) ' s. \color{red} perfect value']);
V_all(1,:)=axis;
hold on
% Plot D1
subplot(2,1,2);
hist11=histogram(D1,'normalization','probability','BinLimits',BinLimits,'BinWidth',BinWidth,'NumBins',NumBins);
title(['D only symbol 1 present. tau = ' num2str(tau_now) ' s. \color{red} perfect value']);
V_all(2,:)=axis;
hold on


% Make axis limits identical
V=[min(V_all(:,1)) max(V_all(:,2)) min(V_all(:,3)) max(V_all(:,4))];
axis(V);
subplot(2,1,1);
axis(V);

% Plot perfect results
subplot(2,1,1);
plot(repmat(D0_perfect,1,2),V_use(3:4),'r-');
subplot(2,1,2);
plot(repmat(D1_perfect,1,2),V_use(3:4),'r-');

% compute overlapping area, where each distribution area normalized to 0.5 as  per paper 1
overlap_area=0;
for i=1:NumBins
 overlap_area = overlap_area + min([hist00.Values(i) hist11.Values(i)])*0.5;
end

% Paper 1 says BER rate is overlap_area
ber_from_pdf=overlap_area;
subplot(2,1,2);
xlabel(['BER from Stenner method = ' num2str(ber_from_pdf)]);





%-----------------------------------------------------------------------------------------------------------------------
%
%-----------------------------------------------------------------------------------------------------------------------
% Unsure if I will  ever want these plots
if(1==2)

    % Allocate memory for axis limits
    V=zeros(2,4);
    figure(51);
    clf
    % Plot symb_0_in_symbol_0_replica
    ax1=subplot(2,1,1);
    hist(ax1,symb_0_in_symbol_0_replica)
    V(1,:)=axis;
    hold on
    % Plot symb_1_in_symbol_1_replica
    ax2=subplot(2,1,2);
    hist(ax2,symb_1_in_symbol_1_replica)
    V(2,:)=axis;
    hold on
    % Get limits to use for both subplots: V_use
    V_use=[min(V(:,1)) max(V(:,2)) min(V(:,3)) max(V(:,4))];
    % Make axis limits same and draw time direct path99
    for i=1:29
        subplot(2,1,i);
        axis(V_use);
    end

    % Plot value of matched filter when no noise
    subplot(2,1,1);
    plot(repmat(mf_perfect_symb_0_in_symbol_0_replica,1,2),V_use(3:4),'r');
    xlabel("Match Filter Output Symbol 0 replica and symbol 0 input");
    title('\color{red} Matched filt output no noise symbol 0');


    subplot(2,1,2);
    plot(repmat(mf_perfect_symb_1_in_symbol_1_replica,1,2),V_use(3:4),'r');
    xlabel("Match Filter Output Symbol 1 replica and symbol 1 input");
    title('\color{red} Matched filt output no noise symbol 1');
end % if(1==2)