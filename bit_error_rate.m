% Derives the bit-error rate as a function of time following the procedure outlined in Fig. 3 of
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
% n_symbol_transmit    1 x 1. Defined by simulate_c_information_measurement.m
% n_samps_before_c
%                      1 x 1. Calculations of the bit error rate are tabulated starting at
%                      n_samps_before_c samples before the time the point of non-analyticity arrives at receiver
%                      in the absence of temporal interference, finite bandwidth filtering, & and noise.
%                      Is >=0.  Zero corresponds to time when light arrrives at receiver in vacuum when there is
%                      no interference (and not infinite bandwidth and no noise).  Want to make  n_samps_before_c
%                      large enough so you can see if detection of symbols can be discerened before light would get there
%                      in a vacuum. Sample interval is 1/f_samp (above).
% d_rec                Structured array defined by output of  simulate_c_inforrmation_measurement_part1.m which calls it
%                      h.d_rec
% fign                 1 x 1.  Matlab figure for plotting ber_info. If emtpy, no plot made
%
% OUTPUTS
%   
% ber_info             Structured array with fields
%
%  -tau                1 x q.  
%                      Integration time (s) with respect to ind_na_arrive_no_filters_no_interference.
%                      tau(n_samps_before_c + 1) = 0.  Interval is dt.
%                      tau(n_samps_before_c+1) = 0 corresponds to index ind_na_arrive_no_filters_no_interference.
%  -ber                1 x q.  ber(j) is bit error rate at tau(j).
% ierr                 0: no error detected. 1 otherwise
%
function [ber_info,ierr]=bit_error_rate(symbol_0_no_noise,symbol_1_no_noise,dt,t_s,ind_direct_na_arrive,...
    ind_tau_alpha,ind_na_arrive_no_filters_no_interference,snr,n_symbol_transmit,n_samps_before_c,d_rec,fign)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preliminary Steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
progname='bit_error_rate.m';
% fprintf(1,'%s: See jls7ub.mat to debug. pausing\n',progname);save jls7ub;pause


% Initialize outputs
ber_info.tau=[];
ber_info.ber=[];
ierr=0;
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

% Derive factor for output of filt_the_data to produce standard deviation of std_n: factor_to_get_std_n
[xx]=filt_the_data(d_rec,randn(1,4*nt),4*nt);
factor_to_get_std_n = std_n/std(xx);

% Initialize
alpha_0_tau_alpha=0;
alpha_1_tau_alpha=0;
% Derive alpha_0_tau_alpha & alpha_1_tau_alpha by using mean values over n_symbol_transmit realizations of noise
for i=1:n_symbol_transmit
    % Generate noise for both symbols
    n_for_symbs=filt_the_data(d_rec,randn(2,nn),[nn nn]);
    n_for_symbs=factor_to_get_std_n * n_for_symbs;

    % Generate noisy data for symbol_0_no_noise & symbol_1_no_noise: chi)0 & chi_1
    chi_0=symbol_0_no_noise(ind_t_s:ind_tau_alpha) + n_for_symbs(1,:);
    chi_1=symbol_1_no_noise(ind_t_s:ind_tau_alpha) + n_for_symbs(2,:);
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

% For debugging, set index of second dimension of D0 & D1 to store values for output of matched filter
ind_D_store=25;
symb_0_in_symbol_0_replica=nan(1,n_symbol_transmit);
symb_0_in_symbol_1_replica=nan(1,n_symbol_transmit);
symb_1_in_symbol_0_replica=nan(1,n_symbol_transmit);
symb_1_in_symbol_1_replica=nan(1,n_symbol_transmit);

for i=1:n_symbol_transmit
    if(rem(i,10)==0)
        fprintf(1,'Starting simulation %d of %d\n',i,n_symbol_transmit);
    end
    % Generate noise for both symbols
    n_for_symbs=filt_the_data(d_rec,randn(2,nt),[nt nt]);
    n_for_symbs=factor_to_get_std_n * n_for_symbs;
    % Generate noisy data for symbol_0_no_noise & symbol_1_no_noise: chi_0 & chi_1
    chi_0=symbol_0_no_noise + n_for_symbs(1,:);
    chi_1=symbol_1_no_noise + n_for_symbs(2,:);
    % chi_0=symbol_0_no_noise;
    % chi_1=symbol_1_no_noise;

    
    % Loop over ind_now
    % for ind_now=ind_tau_alpha:nt
    for ind_now=ind_first:nt
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
       
        % Store data for ind_D_store
        if(ind_D==ind_D_store)
            symb_0_in_symbol_0_replica(i)=l_0;
            symb_0_in_symbol_1_replica(i)=l_1;
        end
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
        % Store data for ind_D_store
        if(ind_D==ind_D_store)
            symb_1_in_symbol_0_replica(i)=l_0;
            symb_1_in_symbol_1_replica(i)=l_1;
        end
    end % for ind_now=ind_tau_alpha:nt

end % for i=1:n_symbol_transmit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Debugging: Plot histogram of noisy symbol 0 output from symbol  0 and symbol 1 matched filter
%            Plot histogram of noisy symbol 1 output from symbol  0 and symbol 1 matched filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
debugplt=0;
if(debugplt==1)
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


    % Get limits to use for both figures: V_use
    V_use=[min(V(:,1)) max(V(:,2)) min(V(:,3)) max(V(:,4))];
    % make both figure have same axis limits
    for i=1:2
        subplot(2,1,i);
        axis(V_use);
    end
    figure(50);
    for i=1:2
        subplot(2,1,i);
        axis(V_use);
    end
end % if(debugplt=1)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Debugging: Plot D0(:,j) and D1(:,j) as functions of their index, j
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
debugplt=0;
if(debugplt==1)
    % Allocate memory to store axes
    V_all=zeros(2,4);
    % index when earliest point of non-analyticity occurs for direct path. earliest is over both symbols
    d_ind = ind_direct_na_arrive - ind_tau_alpha;
    % When tau below occurs for j=1, what is
    for j=1:nD
        % Compute value of tau corresponding to j
        % NO longer true:  tau=0 occurs at index ind_tau_alpha
        % Now true: tau=0 corresponds to index ind_na_arrive_no_filters_no_interference
        % tau_now=(j-1)*dt;
        tau_now=( (j - 1) - n_samps_before_c)*dt;

        figure(100);
        clf;
        subplot(2,1,1);
        hist0=histogram(D0(:,j),'normalization','probability');
        xlabel('\tau (s)');
        title(['D only symbol 0 present. tau = ' num2str(tau_now) ' s']);
        V_all(1,:)=axis;
        subplot(2,1,2);
        hist1=histogram(D1(:,j),'normalization','probability');
      
        title(['D only symbol 1 present. tau = ' num2str(tau_now) ' s']);

        V_all(2,:)=axis;

        % Make axis limits identical
        V=[min(V_all(:,1)) max(V_all(:,2)) min(V_all(:,3)) max(V_all(:,4))];
        axis(V);
        subplot(2,1,1);
        axis(V);



    end % for j=1:nD
end % if(debugplt==1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit Gaussian distribution of D0(:,j) & D1(:,j) for each j
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allocate memory to store mean and standard deviation of D0(:,j) & D1(:,j) for each j
% D0_mean_std(j,1)  is mean               of D0(:,j)
% "             2   is standard deviation of "
% D1_mean_std is same except pertains to D1
D0_mean_std=nan(nD,2);
D1_mean_std=nan(nD,2);

% Fill them
for j=1:nD
    D0_mean_std(j,:)=[mean(D0(:,j)) std(D0(:,j))];
    D1_mean_std(j,:)=[mean(D1(:,j)) std(D1(:,j))];
end % for j=1:nD


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derive bit error rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Probality density function is
%
%  f(x) = exp( -(x-x_mean)^2 / (2*x_std^2) ) / (sqrt(2*pi*x_std^2))     (1)
%
%  and its area is 1.    I verified this numerically. 
%
% Therefore Gaussian distributions with area 0.5 are
%
% exp( -(x-x_mean)^2 / (2*x_std^2) ) / (sqrt(2*pi*x_std^2)) / 2         (2)
%

factor=0.5;
fign_overlap=[];
% Allocate memory to store error from overlapping_areas_two_normal_dist
ierr_all=zeros(1,nD);
% Allocate memory for tau and ber
tau=nan(1,nD);
ber=nan(1,nD);

% This part of code may take the longest.
fprintf(1,'%s: Start computing overlapping_areas_two_normal_dist.m\n',progname);
t_start=tic;
% Fill tau & ber
parfor j=1:nD
    % tau=0 occurs at index ind_tau_alpha.
    % tau(j)=(j-1)*dt;
    tau(j)=( (j - 1) - n_samps_before_c)*dt;
    [ber(j),ierr_all(j)]=overlapping_areas_two_normal_dist(D0_mean_std(j,:),D1_mean_std(j,:),factor,...
        fign_overlap);
end % for j=1:nD
t_elapsed=toc(t_start);
fprintf(1,'%s: Took %f min.\n\n',progname,t_elapsed/60);


% Check for errors
if(sum(abs(ierr_all))~=0)
    fprintf(1,'%s: Error occured for at least one instance of overlapping_areas_two_normal_dist.m\n',progname);
    ierr=1;
    return
end

% Make output
ber_info.tau=tau;
ber_info.ber=ber;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot if wanted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(~isempty(fign))
   figure(fign);
   clf;
   plot(ber_info.tau,ber_info.ber,'k-');
   xlabel('\tau (s)');
   ylabel('BER');
end


