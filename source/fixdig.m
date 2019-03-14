function y = fixdig(x, n)
%FIXDIG Round towards zero with a specified number of digits.
%
%   Y = FIXDIG(X, N) rounds the elements of X to N digits.
%
%   For instance, fixdig(10*sqrt(2) + i*pi/10, 4) returns 14.14 + 0.3141i
%
%   See also: FIX, FLOOR, CEIL, ROUND, FIXDEC, ROUNDDIG, ROUNDDEC.

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


%   Author:      Peter J. Acklam
%   Time-stamp:  2004-09-22 20:07:59 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam
% ---
% $Id$
% $Date: 2014-01-31 15:27:57 -0500 (Fri, 31 Jan 2014) $
% $LastChangedBy: schum $
% $Rev$

   % Check number of input arguments.
   error(nargchk(2, 2, nargin));

   % Quick exit if either argument is empty.
   if isempty(x) || isempty(n)
      y = [];
      return
   end

   % Get size of input arguments.
   size_x   = size(x);
   size_n   = size(n);
   scalar_x = all(size_x == 1);         % True if x is a scalar.
   scalar_n = all(size_n == 1);         % True if n is a scalar.

   % Check size of input arguments and assign output argument.
   if ~scalar_x && ~scalar_n && ~isequal(size_x, size_n)
      error([ 'When both arguments are matrices they must have' ...
               ' the same size' ]);
   end

   % Real part of X.
   k = find(real(x));
   if ~isempty(k)
      xreal = real(x(k));
      m     = nextpowof10(xreal);
      if scalar_x                       % X is scalar.
         f = 10.^(n - m);
         y = fix(xreal .* f) ./ f;
      else
         y = zeros(size_x);
         if scalar_n                    % N is scalar, X is not.
            f = 10.^(n - m);
         else                           % Neither X nor N is scalar.
            f = 10.^(n(k) - m);
         end
         y(k) = fix(xreal .* f) ./ f;
      end
   end

   % Imaginary part of X.
   k = find(imag(x));
   if ~isempty(k)
      ximag = imag(x(k));
      m = nextpowof10(ximag);
      if scalar_x                       % X is scalar.
         f = 10.^(n - m);
         y = y + i*fix(ximag .* f) ./ f;
      else
         if scalar_n                    % N is scalar, X is not.
            f = 10.^(n - m);
         else                           % Neither X nor N is scalar.
            f = 10.^(n(k) - m);
         end
         y(k) = y(k) + i*fix(ximag .* f) ./ f;
      end
   end
