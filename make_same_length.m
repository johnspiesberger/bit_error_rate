% Given two vectors, outputs same two vectors except output of each is same length equaling 
% longer of two input vectors.  Zeros are padded at end of shorter input.
%
% v1            1 x n vector
% v2            1 x m vector
%
% OUPUTS
%
% x1            1 x q: x1(1:n) = v1 with zeros padded at end if needed. q=max([n m])
% x2            1 x q: x2(1:m) = v2 "
%
function [x1,x2]=make_same_length(v1,v2)


%-------------------------------------------------------------------------------------------------------------------
% Make v1 and v2 same length
%-------------------------------------------------------------------------------------------------------------------
n=length(v1);
m=length(v2);
n_max=max([n m]);

if(n==m)
    x1=v1;
    x2=v2;
else
    
    if(n<m)
        % Zero padd v1
        x1=[v1 zeros(1,n_max-n)];
        x2=v2;
    else
        % n > m
        x1=v1;
        % Zero pad v2
        x2=[v2 zeros(1,n_max-m)];
        
    end
end % if(n==m)