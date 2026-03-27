% Formats a c3d figure and two figures from plot_loci_pair_cst_c.m into a figure with 3x1 subplots
%     - subplot 1: c3d field from c3d_exact_many_cases
%     - subplot 2: loci with tdoa calculated via peak to peak method
%     - subplot 3:      "       "       "        cross-correlation    
%
%
% INPUTS
%
% fnam_c3d  1 x 1.  String file name of a .fig file for the c3d figure in the vertical plane (no extension)
%
% fnam_p2p  1 x 1.   "       "       "       "       "      loci     "       "       "       "       "      with 
%                   tdoa calcuated as the difference in arrical times of the first energetic peaks.
%
% fnam_ccf  1 x 1.  String file name of a .fig file for the loci figure in the vertical plane (no extension) with 
%                   tdoa calculated as the peak time of the ccf between the two signals
%
% fnam_out  1 x 1.  String file name which  the outputted .fig and .jpeg files will be given (sans extension).
%
% fign      1 x 1.  Figure number for plot. May not be empty.
%
% OUTPUTS
%
% ierr      0: No error detected. 1: otherwise

function [ierr] = make_xz_c3d_isd_hyp_3by1(fnam_c3d,fnam_p2p,fnam_ccf,fnam_out,fign)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRELIMINARY STEPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
progname='make_xz_c3d_isd_hyp_3by1.m';
%fprintf (1, '%s: See cags7031.mat to debug.  Pausing\n', progname);save cags7031; pause

%initialize outputs
ierr = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare the child figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------------------------------------------------------------------------------------------------------
%open and extract the input figures via their filenames
%-----------------------------------------------------------------------------------------------------------------------
%make a cell array of the filenames for extraction and create handles
fnams = {fnam_c3d; fnam_p2p; fnam_ccf};
axs = gobjects(3,1);

for k = 1:3
    fig = openfig(fnams{k}, 'invisible');
    ax = findall(fig, 'Type', 'axes');
    axs(k) = ax;
end
%-----------------------------------------------------------------------------------------------------------------------
% plot the receivers on the c3d figure
%-----------------------------------------------------------------------------------------------------------------------
%make the c3d axis current and find c3d value needed to be above surface
axes(axs(1));
hold on;
h_surf = findobj(gca, 'Type', 'Surface');
c3d_data = get(h_surf, 'ZData');
zval = max(c3d_data(:)) + 1;
%--- 2. Loop over the two source axes and copy markers -------------------
for srcAx = axs(2:3)
    % find every hollow-circle marker (each is a single point)
    hLines = findall(srcAx, 'Type','line', 'Marker','o');
    for h = hLines.'
        % each line has one point, so grab the first X and Y
        x = h.XData(1);
        y = h.YData(1);
        % plot that point on the surface axes at height zval
        plot3(axs(1), x, y, zval, 'ko', 'MarkerSize',10,'LineWidth',1);
    end
end

%--- 3. Keep the view flattened (XZ plane) if desired --------------------
view(axs(1), 2);    % comment this out if you want to stay in 3-D
hold off;
%-----------------------------------------------------------------------------------------------------------------------
% apply the proper lables and limits to the child axes
%-----------------------------------------------------------------------------------------------------------------------
%the x and z bounds should already be the same
xlim(axs(1), 'tight');
x_bnds = xlim(axs(1));
ylim(axs(1), 'tight');
z_bnds = ylim(axs(1));
z_bnds(2) = 0;
for i = 1:3
    ax = axs(i);
    axes(ax);
    xlabel('')
    zlabel('z (m)');
    if i == 3
        xlabel('x (m)')
    end
    xlim(ax, x_bnds)
    ylim(ax, z_bnds)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make the parent figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%create the figure
figure(fign);
clf;
%make the subplot position vectors
pos = NaN(3,4);
%manually set x, width, and height
pos(:,1) = 0.11;
pos(:,3) = 0.68;
pos(:,4) = 0.24;
%iterate through 3x1 labeled subplots to get the y
for i = 1:3
    h_sp=subplot(3,1,i);
    title(h_sp,'');
    xlabel(h_sp,'x (m)')
    ylabel(h_sp,'z (m)')
    colorbar(h_sp);
    pos_i = get(gca,'Position');
    pos(i,2) = pos_i(2)+ 0.012*(4-i);
    delete(h_sp)
end
%initialize a list for the new axes
axs_new = gobjects(3,1);
%add the subplots to the parent (with text identifiers A, B, C)
letters = [{'A'},{'B'},{'C'}];
for i = 1:3
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
    %label the subplots
    letter = letters(i);
    letter_x = 0.95*(x_bnds(2) - x_bnds(1))+x_bnds(1);
    letter_z = 0.9*(z_bnds(2) - z_bnds(1))+z_bnds(1);
    text(letter_x,letter_z,letter,'FontSize',15,'FontWeight','bold','Parent',new_ax);
    %save handle
    axs_new(i) = new_ax;
end
% Link x and y axes across all subplots, then save the parent figure
linkaxes(axs_new, 'xy');
saveas(figure(fign),fnam_out,'fig');
exportgraphics(figure(fign),[fnam_out '.jpeg'],'Resolution',600);
%close the temporary figures
for i = 1:3
    ax = axs(i);
    axes(ax);
    close(gcf);
end

