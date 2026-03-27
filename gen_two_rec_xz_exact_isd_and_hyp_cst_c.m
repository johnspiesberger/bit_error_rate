% Runs c3d_exact_many_cases.m and plot_loci_rec_pair_cst_c.m to make a figure with 5 x 1 subplots:
%     - subplot 1: c3d field from c3d_exact_many_cases
%     - subplot 2: hyperbola ~51m from source at closest point, tdoa = 0.0913 s calculated via peak to peak method
%     - subplot 3: no hyperbola exists,                         "      0.1212   "       "       "       "       "       
%     - subplot 4: hyperbola ~66m from source at closest point, "      0.0866   "       "      cross correlation
%     - subplot 5: no hyperbola exists,                         "      0.1251   "       "       "       "             
% It then runs make_xz_c3d_isd_hyp_3by1.m to make 2 figures (A and B) with 3 x 1 subplots:
%     - subplot 1: same as subplot 1 above
%     - subplot 2:  "       "      2    "  for A; same as subplot 4 above for B
%     - subplot 3:  "       "      3    "       "       "       " 5 "       " 
%
% REPOSITORIES USED
%
% effective_speed_interference/dispersion_relation/mlab_simulations         changeset 3
%
% NO INPUTS
%
% OUTPUTS
% out_st    Structured array with fields:
%
%  - fnams  Structured array with fields:
%
%     - c3d       1 x 1.  String file name (sans extension) of a .fig file for the figure of c3d in the vertical plane 
%     - p2p_A     1 x 1.  "       "       "       "       "       "       "        loci figure    "       "       " 
%                         with tdoa calcuated as the difference in arrical times of the first energetic peaks of the 
%                         signal arriving at the receivers in configuration A (see 'hardcode primary input parameters' 
%                         'preliminary steps)
%     - p2p_B     1 x 1.  same as fnam_p2p_A except for receivers in configuration B
%
%     - ccf_A     1 x 1.  String file name (sans extension) of a .fig file for the loci figure in the vertical plane  
%                         with tdoa calculated as the peak time of the ccf between the signals arriving at the receivers 
%                         in configuration A (see 'hardcode primary input parameters' in 'preliminary steps)
%     - ccf_B     1 x 1.  same as fnam_ccf_A except for receivers in configuration B
%     - all       1 x 1.  String filename (sans extension) of a .fig and .jpeg file for the main 5-by-1 subplot
%
%  - data   Structured array with fields:
%
%     - c3d_lims  1 x 2. The min (c3d_lims(1)) and max (c3d_lims(2)) of the c3d (m/s) field
%     - p2p       2 x 2. For receivers in configuration i, p2p(i,1) is the tdoa using the peak to peak method and
%                        p2p(i, 2) is the minimum distance between the source and the hyperbola at the closest point.
%                        p2p(i, 2) = NaN means no hyperbola exists for a given tdoa and receiver configuration
%     - ccf       2 x 2. Same as p2p but for tdoa calculated via ccf
%
% ierr      0: No error detected. 1: otherwise

function [out_st, ierr] = gen_two_rec_xz_exact_isd_and_hyp_cst_c()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRELIMINARY STEPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
progname='gen_two_rec_xz_exact_isd_and_hyp_cst_c.m';
%fprintf (1, '%s: See cags6231.mat to debug.  Pausing\n', progname);save cags6131; pause

%initialize outputs
ierr = 0;
%-----------------------------------------------------------------------------------------------------------------------
% hardcode primary input parameters
%-----------------------------------------------------------------------------------------------------------------------
% set source and receiver coordinates. Configuration "A" is subplots 2 and 4, and "B" is for subplots 3 and 6
src_xz=[5 -50];
rec_1_xz_A=[100 -70];
rec_2_xz_A=[250 -80];
rec_1_xz_B=[80 -70];
rec_2_xz_B=[250 -70];

%set signal parameters
tau=0.025;
omega=93*2*pi;
g_base=2000;
g_mult = 30;
cheby = 5;
mult_type = 0;
boundary_type=1;

%set remaining parameters
fign = 14;
c = 1500;
x_bnds = [-75 275];
z_bnds = [-100 -1];
xz_res = 150;
td_type_p2p = 0;
td_type_ccf = 1;
fign_empty = [];

%set a unique id to append to the figures' filenames to differentiate between runs with the same xz_res and g values
unique_id = '';

%-----------------------------------------------------------------------------------------------------------------------
% calculate secondary input parameters
%-----------------------------------------------------------------------------------------------------------------------
% derive sample frequency
g = g_mult * g_base;
fs = g*omega/(2*pi);

% derive td for configuration A
[h1A,~]=exact_time_domain_solution_reflected_direct_paths(c,src_xz,rec_1_xz_A,tau,omega,g,boundary_type,fign_empty);
[h2A,~]=exact_time_domain_solution_reflected_direct_paths(c,src_xz,rec_2_xz_A,tau,omega,g,boundary_type,fign_empty);
td_p2p_A = h2A.t_r - h1A.t_r;
[td_ccf_A,~]=plot_ah_ccf_same_samp_freq(h1A.s,h2A.s,fs,fign_empty);


% derive td for configuration B
[h1B,~]=exact_time_domain_solution_reflected_direct_paths(c,src_xz,rec_1_xz_B,tau,omega,g,boundary_type,fign_empty);
[h2B,~]=exact_time_domain_solution_reflected_direct_paths(c,src_xz,rec_2_xz_B,tau,omega,g,boundary_type,fign_empty);
td_p2p_B = h2B.t_r - h1B.t_r;
[td_ccf_B,~]=plot_ah_ccf_same_samp_freq(h1B.s,h2B.s,fs,fign_empty);

%-----------------------------------------------------------------------------------------------------------------------
% initialize outputs
%-----------------------------------------------------------------------------------------------------------------------
out_st = struct();

% Substructure: fnams (all are strings, initialized to empty)
out_st.fnams = struct( ...
    'c3d',   '', ...
    'p2p_A', '', ...
    'p2p_B', '', ...
    'ccf_A', '', ...
    'ccf_B', '', ...
    'all',   '' ...
);

% Substructure: data
out_st.data = struct( ...
    'c3d_lims', [NaN, NaN], ...  
    'p2p', NaN(2,2), ...      
    'ccf', NaN(2,2) ...          
);

ierr = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make the individual subplot figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the filename extesion tag
if ~isempty(unique_id)
    unique_id = ['_' unique_id];
end
fnam_tag = ['_res' num2str(xz_res) '_g' num2str(g_base) 'm' num2str(g_mult) unique_id];

%-----------------------------------------------------------------------------------------------------------------------
%make the c3d figure
%-----------------------------------------------------------------------------------------------------------------------
fprintf(1,'%s: Making c3d figure\n',progname)
%set parameters then generate the figure
fign_1 = 624511;
rec_x_bnds = x_bnds;
rec_z_bnds = z_bnds;
nx = xz_res;
nz = xz_res;
c3d_exact_many_cases(c, src_xz, rec_x_bnds, rec_z_bnds, nx, nz, tau, omega, g_base, g_mult, cheby, mult_type, ...
    boundary_type, fign_1);
hold on;
%plot the source as filled circle above the z-level of the surface max
h_surf = findobj(gca, 'Type', 'Surface');
c3d_data = get(h_surf, 'ZData');
c3d_max = max(c3d_data(:));
c3d_min = min(c3d_data(:));
cval = c3d_max + 1e-3*(c3d_max - c3d_min);
plot3(src_xz(1), src_xz(2), cval, 'k.', 'MarkerSize', 27);
%change the view to the 2d xz plane
view(gca, 2);
%reset the titles
title(gca,'');
xlabel(gca, '')
ylabel(gca, 'z (m)')
%remove the source line and minimum marker
delete(findall(fign_1, 'Type', 'line', 'Marker', '*', 'Color', 'k'));  
delete(findall(fign_1, 'Type', 'line', 'LineStyle', '-', 'Color', 'k'));
hold off;
ax1 = gca;
saveas(figure(fign_1),['xz_c3d' fnam_tag],'fig');

%-----------------------------------------------------------------------------------------------------------------------
%make the first peak-to-peak loci figure
%-----------------------------------------------------------------------------------------------------------------------
fprintf(1,'%s: Making peak-to-peak loci figure 1 of 2\n',progname)
%generate the figure
fign_2 = 624512;
[min_dist_p2p_A,~] = plot_loci_rec_pair_cst_c(c, src_xz, rec_1_xz_A, rec_2_xz_A, tau, omega, g_base, g_mult, cheby, ...
    mult_type, boundary_type,  x_bnds, z_bnds, xz_res, td_p2p_A, td_type_p2p, fign_2);
hold on;
%reset the titles
title(gca,'');
xlabel(gca, '')
ylabel(gca, 'z (m)')
%delete the legend and any previously plotted sources or receivers
lgd = legend(gca); 
delete(lgd);  
delete(findall(fign_2, 'Type', 'line', 'Marker', 'o'));
%plot the source as filled circle and receivers as hollow circles
plot(src_xz(1),src_xz(2),'k.','MarkerSize',27)
plot(rec_1_xz_A(1),rec_1_xz_A(2),'ko','MarkerSize',10,'LineWidth',1)
plot(rec_2_xz_A(1),rec_2_xz_A(2),'ko','MarkerSize',10,'LineWidth',1)
hold off;
%store the axis handle
ax2 = gca;

saveas(figure(fign_2),['xz_isd_hyp_pk2pk_A' fnam_tag],'fig');

%-----------------------------------------------------------------------------------------------------------------------
%make second peak-to-peak loci figure
%-----------------------------------------------------------------------------------------------------------------------
fprintf(1,'%s: Making peak-to-peak loci figure 2 of 2\n',progname)
%generate the figure
fign_3 = 624513;
[min_dist_p2p_B,~] = plot_loci_rec_pair_cst_c(c, src_xz, rec_1_xz_B, rec_2_xz_B, tau, omega, g_base, g_mult, cheby, ...
    mult_type, boundary_type, x_bnds, z_bnds, xz_res, td_p2p_B, td_type_p2p, fign_3);
hold on;
%reset the titles
title(gca,'');
xlabel(gca, 'x (m)')
ylabel(gca, 'z (m)')
%delete the legend and any previously plotted sources or receivers
lgd = legend(gca); 
delete(lgd);  
delete(findall(fign_3, 'Type', 'line', 'Marker', 'o'));
%plot the source as filled circle and receivers as hollow circles
plot(src_xz(1),src_xz(2),'k.','MarkerSize',27)
plot(rec_1_xz_B(1),rec_1_xz_B(2),'ko','MarkerSize',10,'LineWidth',1)
plot(rec_2_xz_B(1),rec_2_xz_B(2),'ko','MarkerSize',10,'LineWidth',1)
hold off;
ax3 = gca;
saveas(figure(fign_3),['xz_isd_hyp_pk2pk_B' fnam_tag],'fig');

%-----------------------------------------------------------------------------------------------------------------------
% make first ccf loci figure
%-----------------------------------------------------------------------------------------------------------------------
fprintf(1,'%s: Making ccf loci figure 1 of 2\n',progname)
%generate the figure
fign_4 = 624514;
[min_dist_ccf_A,~] = plot_loci_rec_pair_cst_c(c, src_xz, rec_1_xz_A, rec_2_xz_A, tau, omega, g_base, g_mult, cheby, ...
    mult_type, boundary_type, x_bnds, z_bnds, xz_res, td_ccf_A, td_type_ccf, fign_4);
hold on;
%reset the titles
title(gca,'');
xlabel(gca, '')
ylabel(gca, 'z (m)')
%delete the legend and any previously plotted sources or receivers
lgd = legend(gca); 
delete(lgd);  
delete(findall(fign_4, 'Type', 'line', 'Marker', 'o'));
%plot the source as filled circle and receivers as hollow circles
plot(src_xz(1),src_xz(2),'k.','MarkerSize',27)
plot(rec_1_xz_A(1),rec_1_xz_A(2),'ko','MarkerSize',10,'LineWidth',1)
plot(rec_2_xz_A(1),rec_2_xz_A(2),'ko','MarkerSize',10,'LineWidth',1)
hold off;
ax4 = gca;
saveas(figure(fign_4),['xz_isd_hyp_ccf_A' fnam_tag],'fig');

%-----------------------------------------------------------------------------------------------------------------------
% make second ccf loci figure
%-----------------------------------------------------------------------------------------------------------------------
fprintf(1,'%s: Making ccf loci figure 2 of 2\n',progname)
%generate the figure
fign_5 = 624515;
[min_dist_ccf_B,~] = plot_loci_rec_pair_cst_c(c, src_xz, rec_1_xz_B, rec_2_xz_B, tau, omega, g_base, g_mult, cheby, ...
    mult_type, boundary_type, x_bnds, z_bnds, xz_res, td_ccf_B, td_type_ccf, fign_5);
hold on;
%reset the titles
title(gca,'');
xlabel(gca, 'x(m)')
ylabel(gca, 'z (m)')
%delete the legend and any previously plotted sources or receivers
lgd = legend(gca); 
delete(lgd);  
delete(findall(fign_5, 'Type', 'line', 'Marker', 'o'));
%plot the source as filled circle and receivers as hollow circles
plot(src_xz(1),src_xz(2),'k.','MarkerSize',27)
plot(rec_1_xz_B(1),rec_1_xz_B(2),'ko','MarkerSize',10,'LineWidth',1)
plot(rec_2_xz_B(1),rec_2_xz_B(2),'ko','MarkerSize',10,'LineWidth',1)
hold off;
ax5 = gca;
saveas(figure(fign_5),['xz_isd_hyp_ccf_B' fnam_tag],'fig');

%combine the axes into a vector
axs = [ax1 ax2 ax3 ax4 ax5];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make the parent figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fill the parent figure with the child subplots
figure(fign); clf;
set(gcf, 'Units', 'pixels');
fig_pos = get(gcf, 'Position');
fig_pos(4) = 550; % ← increase height to allow for no overlaps
set(gcf, 'Position', fig_pos);
axs_new = gobjects(1,5);

%make the subplot position vectors
pos = NaN(5,4);
%manually set x, width, and height
pos(:,1) = 0.11;
pos(:,3) = 0.68;
pos(:,4) = 0.135;
%iterate through 3x1 labeled subplots to get the y
for i = 1:5
    h_sp=subplot(5,1,i);
    title(h_sp,'');
    xlabel(h_sp,'x (m)')
    ylabel(h_sp,'z (m)')
    colorbar(h_sp);
    pos_i = get(gca,'Position');
    pos(i,2) = pos_i(2) + 0.02*(4-i);
    delete(h_sp)
end

%add the subplots to the parent (with text identifiers A-E))
letters = [{'A'},{'B'},{'C'},{'D'},{'E'}];
for i = 1:5
    % Copy original axis and paste in this figure
    new_ax = copyobj(axs(i), gcf);
    % Adjust for top plot (e.g., c3d field)
    if i == 1
        colormap('gray');
        cb = colorbar(new_ax);
        ylabel(cb, 'c_{3d} (m/s)')
    end
    %set up new axis
    set(new_ax, 'Position', pos(i,:), ...
                'Box', 'on', ...
                'DataAspectRatioMode', 'auto', ...
                'PlotBoxAspectRatioMode', 'auto', ...
                'Layer', 'top', ...
                'FontSize', 12);  % ensure ticks on top of graphics

    % Enforce x and z limits
    xlim(new_ax, x_bnds);
    ylim(new_ax, z_bnds);
    %label the subplots
    letter = letters(i);
    letter_x = 0.95*(x_bnds(2) - x_bnds(1))+x_bnds(1);
    letter_z = 0.8*(z_bnds(2) - z_bnds(1))+z_bnds(1);
    text(letter_x,letter_z,letter,'FontSize',15,'FontWeight','bold','Parent',new_ax);
    %save handle
    axs_new(i) = new_ax;
end

% Link x and y axes across all subplots, then save the parent figure
linkaxes(axs_new, 'xy');
saveas(figure(fign),['xz_c3d_isd_hyp_cst_c' fnam_tag],'fig');
exportgraphics(figure(fign),['xz_c3d_isd_hyp_cst_c' fnam_tag '.jpeg'],'Resolution',600);
close([fign_1, fign_2, fign_3, fign_4, fign_5])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fill out_st
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Substructure: fnams (all are strings, initialized to empty)
out_st.fnams = struct( ...
    'c3d',   ['xz_c3d' fnam_tag], ...
    'p2p_A', ['xz_isd_hyp_pk2pk_A' fnam_tag], ...
    'p2p_B', ['xz_isd_hyp_pk2pk_B' fnam_tag], ...
    'ccf_A', ['xz_isd_hyp_ccf_A' fnam_tag], ...
    'ccf_B', ['xz_isd_hyp_ccf_B' fnam_tag], ...
    'all',   ['xz_c3d_isd_hyp_cst_c' fnam_tag] ...
    );

% Substructure: data
out_st.data = struct( ...
    'c3d_lims', [c3d_min, c3d_max], ...  % 1x2 vector
    'p2p', [td_p2p_A, min_dist_p2p_A; td_p2p_B, min_dist_p2p_B], ...         % 2x2 matrix
    'ccf', [td_ccf_A, min_dist_ccf_A; td_ccf_B, min_dist_ccf_B] ...          % 2x2 matrix
);
% save the structure data
save(['DAT_xz_c3d_isd_hyp_cst_c' fnam_tag],'out_st')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make the 3 by 1 plot if desired (use plot_3by: 1 x 1. true: plots, false: does not plot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_3by1 = true;

if plot_3by1
    %set filenames
    fnam_c3d = out_st.fnams.c3d;
    fnam_p2p_A = out_st.fnams.p2p_A;
    fnam_p2p_B = out_st.fnams.p2p_B;
    fnam_ccf_A = out_st.fnams.ccf_A;
    fnam_ccf_B = out_st.fnams.ccf_B;
    fnam_out_A = ['FIG_xz_c3d_isd_hyp_A' fnam_tag];
    fnam_out_B = ['FIG_xz_c3d_isd_hyp_B' fnam_tag];
    %set fiure numbers
    fign_A = 7111;
    fign_B = 7112;
    %make 3by figures
    make_xz_c3d_isd_hyp_3by1(fnam_c3d,fnam_p2p_A,fnam_ccf_A,fnam_out_A,fign_A);
    make_xz_c3d_isd_hyp_3by1(fnam_c3d,fnam_p2p_B,fnam_ccf_B,fnam_out_B,fign_B);
end %plot_3by1
