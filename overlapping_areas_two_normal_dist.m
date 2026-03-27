% Computes overlapping area of two normal probability distribution functions
%
% INPUTS
%
% mean_std1        1 x 2.  mean_std1(j) is mean & standard deviation of first normal probability dist function with mean in 
%                  in j=1.
% mean_std2        1 x 2. Same as mean_std1 except pertains to second normal probablity dist function.
% factor           1 x 1. Factor to multiply area of each normal distribution by before deriving overlapping area.
%                  Area is one for each input distribution.
%                  E.g. factor=0.5 means area of each inputted normal distribution is renormalized to 0.5 before
%                  computing overlapping area
% fign             1 x 1. Matlab figure # to show both distributions (times factor) and overlapping area.
%                  If empty, plot not shown.
%
% OUTPUTS
%
% area_overlap     1 x 1. Overlapping area
% ierr             0: no errors detected. 1 otherwise
%
function [area_overlap,ierr]=overlapping_areas_two_normal_dist(mean_std1,mean_std2,factor,fign)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preliminary Steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
progname='overlapping_area_two_normal_dist.m';
% fprintf(1,'%s: See jls7u.mat to debug. pausing\n',progname);pause

dotest=0;
if(dotest==1)
   fprintf(1,'%s: Test case. Press ENTER to continue\n',progname);pause 
   mean_std1=[0 1];
   mean_std2=[.5 2];
   factor=0.5;
   fign=33;
end

% Initialize outputs
area_overlap=NaN;
ierr=0;

% Specify number standard deviations from mean to use for integration
nstd=4;

% Specify fraction of smallest standard deviation to use for integrating on independent axis
f=0.01;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derive independent variable bounds for integration: x_bnds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get limits of from mean_std1
dlim1=mean_std1(1) + nstd*mean_std1(2)*[-1 1];
% Get limits of from mean_std2
dlim2=mean_std2(1) + nstd*mean_std2(2)*[-1 1];
% Get limits of integration
x_bnds=[ min([dlim1(1) dlim2(1)])  max([dlim1(2) dlim2(2)]) ];

% Derive x interval for integration
dx=f*min([mean_std1(2) mean_std2(2)]);

% Derive x axis values to evaluate both functions
x=x_bnds(1):dx:x_bnds(2);

% Get number x: nx
% It might be very large and exceed PC's ram when standard deviations are very small.
nx=length(x);
% Allocate memory to store minimum of both distributions at each x
pdf_min=nan(1,nx);

% Derive first pdf values to use for each x: pdf1
pdf1=factor * exp( -(x-mean_std1(1)).^2 / (2*mean_std1(2)^2) )/ (sqrt(2*pi*mean_std1(2)^2));
% Derive second pdf values to use for each x: pdf2
pdf2=factor * exp( -(x-mean_std2(1)).^2 / (2*mean_std2(2)^2) )/ (sqrt(2*pi*mean_std2(2)^2));

% Compute pdf_min: Its min of pdf1 & pdf2 at each x
pdf_min=min([pdf1;pdf2]);

% Compute overlapping area
area_overlap=trapz(pdf_min)*dx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot if wanted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(~isempty(fign))
   figure(fign);
   clf;
   plot(x,pdf1,'k-',x,pdf2,'g--',x,pdf_min,'r-');
   xlabel('x');
   title(['pdf1 \color{green} pdf2 \color{red} overlap area = ' num2str(area_overlap)]); 
   fprintf(1,'%s: Pausing. Press ENTER to continue\n',progname);
   pause
   close(fign);
end


