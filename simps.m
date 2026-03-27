% simps.m computes the integral of a function using simpson's rule.
% [ans1,err] = simps(y,a,b) where the integrand is in the vector y with an odd
% number of entries, and the integrand falls between the independent 
% variable limits from a to b.  John Spiesberger 17 Oct 1993.
% The integral is output in ans1.  err=0 if ok, and not equal to 0
% if integral has an error.
%     $RCSfile: simps.m,v $
%     $Revision: 1.2 $
%     $Date: 2003/05/07 13:06:05 $
%     $Author: jspies $
  function [ans1,err] = newsimps(y,a,b)
  err = 0;
  [nrow,ncol] = size(y);
  n = max(nrow,ncol);
  if( n < 3 )
    err = 1;
    return
  end,
  if( rem(n,2) == 0  ),
    % there are an even # of integrands
    err = 2;
    return
  end
  
  ans1 = y(1) + y(n);

  % get contribution from even numbered indicies
  kend = n - 1;
 
  sumeven=sum(y(2:2:kend));
  ans1 = ans1 + 4. * sumeven;

  if( n > 3 ),
    sumofodd=sum(y(3:2:kend));
    ans1 = ans1 + 2.*sumofodd;
  end;
  
  ans1 = ( (b-a)/(n-1) ) / 3. * ans1;
  return
