% Derives the 3D effective speed for many cases of exact_time_domain_solution_surf_direct_paths.m
%
% INPUTS
%
% c             1 x 1.  speed of sound in water (m/s)
% src_z         1 x 1.  Source is located at  horizontal coordinate x=0 and  vertical cooordinate, z=src_z (m).
%               z is 0 at ocean's surface and is positive up.  source also located at y=0,  in 3D space.
% rec_x_bnds    1 x 2 Receiver located at src_z between min/max values of x (m) with min in col 1.
% rec_z         1 x 1.  Vertical coordinate of receiver (m). Same coordinate system as for src_z
% q             1 x 1.  Receiver x coordinate will be simulated q times in rec_x_bnds
% tau           1 x 1. emitted pulse is exp(-(t/tau)^2)*cos(omega*t) with tau in seconds and t in seconds
% omega         1 x 1. See tau above. Units are rad/s
% boundary_type 1 x 1. 1:  boundary is like surface of ocean where there are no pressure fluctuations.
%                          There is a direct path and a path reflecting from surface 
%                      2: boundary is hard. Reflected path bounces from hard flat rocks on bottom of ocean.
% fign          1 x 1. Figure number for plot. If empty, no plots shown
%
% OUTPUTS
%
% ierr          0: No error detected. 1: otherwise
%
function [ierr]=c3d_exact_many_cases(c,src_z,rec_x_bnds,rec_z,q,tau,omega,boundary_type,fign)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preliminary Steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
progname='c3d_exact_many_cases.m';

dotest=1;
if(dotest==1)
  fprintf(1,'%s: Test case. Press Enter to continue\n',progname);pause
  c=1500;
  scenario=4;
  switch scenario
      case 1
          % Yields c3d < c. i.e. ~900 m/s
          src_z=-4.9;
          rec_x_bnds=[0.1 4];
          rec_z=-5;
          q=100;
          tau=0.025;
          omega=100*2*pi;
          boundary_type=1;
      case 2
          % Yields c3d between 1160 and 2013 m/s
          src_z=-25;
          rec_x_bnds=[0.1 100];
          rec_z=-25;
          q=100;
          tau=0.025;
          omega=100*2*pi;
          boundary_type=1;
      case 4
          
          src_z=-30;
          rec_x_bnds=[0.1 100];
          rec_z=-25;
          q=100;
          tau=0.025;
          omega=100*2*pi;
          boundary_type=1;
          
  end % switch scenario
  fign=1;
end

% Initialize outputs
ierr=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do work
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allocate memory
c3d=nan(1,q);
% Generate x coords of receiver
rec_x=linspace(rec_x_bnds(1),rec_x_bnds(2),q);
fign_exact=[];
% Fill c3d
for i=1:q
    % Derive xz coord of receiver
    rec_xz=[rec_x(i) rec_z];

    [h,ierr]=exact_time_domain_solution_reflected_direct_paths(c,src_z,rec_xz,tau,omega,boundary_type,fign_exact);
    if(ierr~=0),return,end
    if(~isempty(fign_exact))
        fprintf(1,'%s: Pausing. Press ENTER to plot next case.\n',progname);
        pause
    end
    % Store c3d
    c3d(i)=h.c3d;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot if wanted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(~isempty(fign))
    figure(fign);
    clf;
    plot(rec_x,c3d,'k.');
    hold on;
    % Plot in-situ speed
    V=axis;
    plot(V(1:2),[c c],'g-');
    xlabel('Distance (m) between source and receiver');
    ylabel('c3d (m/s)');
    title(['Min/max c3d = ' num2str(min(c3d)) '  ' num2str(max(c3d)) ' (m/s). \color{green} in-situ speed']);
end
