% This subroutine will generate a hyperbola, points of equal difference in distance from two receivers, in the xz plane
%
% INPUTS
%
% rec_1         1 x 2 The x and z coordinate of the first receiver (x,z) (m). x is the horizontal coordinate, and z is 
%                     the vertical coordinate. z is defined positive up and z = 0 at the ocean's surface
% rec_2         1 x 2 The x and z coordinate of the second receiver (x,z) (m)
% dist_diff     1 x 1 The defined difference in distance between all the points on the hyperbola (m)
% x_bnds        1 x 2 The x limits to evaluate points between. The hyperbola can only be generated on points in this x
%                     range (m)
% z_bnds        1 x 2 The z limits to evaluate points between. The hyperbola can only be generated on points in this z
%                     range (m)
% xz_res        1 x 1 The number of points to divide the x and z bounds into. The higher the number the more points
%                     evaluated at
% fign          1 x 1 If empty contour plot will not show, otherwise contour plot will appear automatically
%
% OUTPUTS
%
% x_coords      1 x n The x coordinates of all the points that fall on the hyperbola (m)
% z_coords      1 x n The z coordinates of all the points that fall on the hyperbola. Goes with x_coords
%                     (x_coords(i),z_coords(i)) (m)

function [x_coords,z_coords] = hyperbola_2d_xz(rec_1, rec_2, dist_diff, x_bnds, z_bnds, x_z_res, fign)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRELIMINARY STEPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
progname='hyperbola_2d_xz.m';
%fprintf (1, '%s: See cags6241.mat to debug.  Pausing\n', progname);save cags6241; pause

% Initialize outputs
x_coords = [];
z_coords = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make hyperbola
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract the receiver coordinates
r1x = rec_1(1);
r1z = rec_1(2);
r2x = rec_2(1);
r2z = rec_2(2);

% Create a mesh to check distances at every point on the mesh
x = linspace(x_bnds(1), x_bnds(2), x_z_res);
z = linspace(z_bnds(1), z_bnds(2), x_z_res);
[X, Z] = meshgrid(x, z);

% Calculate the Euclidean distances from each point on the grid to each receiver
distance1 = sqrt((X - r1x).^2 + (Z - r1z).^2);
distance2 = sqrt((X - r2x).^2 + (Z - r2z).^2);

% Calculate the difference in distances at each point
difference = distance2 - distance1;
% Make a contour at the inputted distance in difference
d_level = [dist_diff dist_diff];
if ~isempty(fign)
    figure(fign);  % Open or reuse figure only if specified
    hold on;
    contour(X, Z, difference, d_level, 'k--', 'ShowText', 'off','LineWidth',2);
    hold off;
end

%extract the data regardless of plotting
M = contourc(x, z, difference, d_level);

% Fix the formatting of the M matrix so that the only points are the points that fall in the contour
index_need_clear = M(1,:)==dist_diff;
M(:,index_need_clear) = [];
% From the contour matrix extract the coordinates of the contour
x_coords = M(1,2:length(M));
z_coords = M(2,2:length(M));


