% Derives (x,z) bounds for locations of source and receiver.
%
% INPUTS
%
% tau             1 x 1.  Same as invars.tau       defined by simulate_c_information_measurements.m
% tau_after       1 x 1.  "             .tau_after "
% passband_freq1  1 x 1.  "             .passband_freq1 "
% passband_freq2  1 x 1.  "                       freq2 "
% which_crit      1 x 1.  1:  use crit_dist for critical distance
%                         2:  use crit_dist_optical_switch for critical distance
%                 where
%                 pulse_resol=min([tau tau_after]);
%                 crit_dist = invars.c*pulse_resol/2;
%                 crit_dist_optical_switch = 1/(passband_freq2-passband_freq1)/2;
% f_mult          1 x 1. multiplication factor, f_mult, for setting  bounds of source & receiver location: f_mult
% c               1 x 1. Same as invers.c defined by simulate_c_information_measurements.m
%             
% OUTPUTS
%
% loc_bnds        2 x 2. loc_bnds(1,:) are min/max values of x coordinate for choosing locations of source and receiver.
%                 Min value in col. 1.  x coordinate of source is always 0. Max value >=0
%                 loc_bnds(2,:) are min/max values of z coordinates of source and receiver with min in col. 1.
%                 z positive up and 0 at boundary.  Units meters.
%                 loc_bnds =[0 f_mult*crit_dist_use; ...
%                           -f_mult*crit_dist_use -1.d-6*crit_dist_use ];
%                 where crit_dist_use is one of crit_dist and crit_dist_optical_switch, determined by which_crit.
% crit_dist_use   1 x 1. Critical distance used (m).  Determined by which_crit.
% ierr            0: no errors detected
%                 1: otherwise
%
function [loc_bnds,crit_dist_use,ierr]=derive_loc_bnds(tau,tau_after,passband_freq1,passband_freq2,which_crit,f_mult,c)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preliminary Steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
progname='derive_loc_bnds.m';
% fprintf(1,'%s: See jlsb0v.mat to debug. pausing\n',progname);save jlsb0v;pause


% Initialize outputs
loc_bnds=nan(2,2);
crit_dist_use=[];
ierr=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do work
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derive approximate pulse resolution of emitted signal, not including bandwidth due to non-analyticity
pulse_resol=min([tau tau_after]);
% Derive critical distance from boundary: crit_dist
crit_dist=c*pulse_resol/2;

% Derive critical distance due to bandwidth of optical switch
crit_dist_optical_switch=c/(passband_freq2-passband_freq1)/2;


% Set which critical distance to use: crit_dist_use
if(which_crit==1)
    crit_dist_use=crit_dist;
else
    crit_dist_use=crit_dist_optical_switch;
end
% Specify bounds in Cartesian (x,z) space for locations of source and receiver: loc_bnds.
% loc_bnds(1,:) are min/max values of x with min in col 1.
% loc_bnds(2,:) are min/max values of z with min in col 1.
% x coordinate of source is always 0
% Units are meters, z positive up and zero at boundary.
loc_bnds=[0 f_mult*crit_dist_use; -f_mult*crit_dist_use -1.d-6*crit_dist_use ];

