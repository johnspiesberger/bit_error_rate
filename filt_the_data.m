% Filters the data using filter design in design_bandpass_filter.m
%
% INPUTS
%
% d                  Output of designfilter.m. Digital filter design
% x                  nrec x n Matrix of data to be filtered with filtfilt.m.
%                    There are nrec time series to filter.
% nx                 1 x nrec.  Valid data for x(i,:) are x(i,1:nx).
%
% OUTPUTS
%
% y                  nrec x n. Matrix of filtered data
%
function [y]=filt_the_data(d,x,nx)
% save jlshh
[nrec,n]=size(x);

% Allocate memory
y=zeros(nrec,n);
for i=1:nrec
    [y(i,1:nx(i))]=filter(d,x(i,1:nx(i)));
end
