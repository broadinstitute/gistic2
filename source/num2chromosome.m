function s=num2chromosome(n)
% NUM2CHROMOSOME convert chromosome index to string representation
%   S = num2chromosome(N)
%   Returns a cell array of strings corresponding the chromosome
%   numbers in N. If N is a scalar, a simple string is returned.
%

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

RGI = RefGeneInfo;
s = getChromosomeText(RGI,n);