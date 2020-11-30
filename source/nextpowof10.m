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
% Copyright (c) 2011-2017 Gad Getz, Rameen Beroukhim, Craig Mermel,
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, Gordon Saksena
% All Rights Reserved.
% (See the accompanying LICENSE file for licensing details.)

   error(nargchk(1, 1, nargin));

   if ~isreal(x)
      error('Input must be real.');
   end

   x = abs(x);
   p = ceil(log(x) / log(10));          % estimate
   k = x <= 10.^(p - 1);
   p(k) = p(k) - 1;                     % correction
