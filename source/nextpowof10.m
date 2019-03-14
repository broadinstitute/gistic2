function p = nextpowof10(x)
%NEXTPOWOF10 Next power of 10.
%
%   P = NEXTPOWOF10(X) returns the smallest integer P such that 10^P >= abs(X).
%
%   Essentially, NEXTPOWOF10(X) is the same as CEIL(LOG(ABS(X)) / LOG(10)), but
%   special care is taken to catch round-off errors.
%
%   See also PREVPOWOF2, NEXTPOW.

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


%   Author:      Peter J. Acklam
%   Time-stamp:  2003-11-17 11:38:03 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam
% ---
% $Id$
% $Date: 2014-01-31 15:28:37 -0500 (Fri, 31 Jan 2014) $
% $LastChangedBy: schum $
% $Rev$

  
  
   error(nargchk(1, 1, nargin));

   if ~isreal(x)
      error('Input must be real.');
   end

   x = abs(x);
   p = ceil(log(x) / log(10));          % estimate
   k = x <= 10.^(p - 1);
   p(k) = p(k) - 1;                     % correction
