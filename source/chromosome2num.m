function n=chromosome2num(ch)
%CHROMOSOME2NUM map chromosome strings to chromosome index
%
%   N=chromosome2num(CH)
%
%   replace a chromosome string or array of strings by numbers
%   invalid strings generate NaNs 

% GISTIC software version 2.0
% Copyright (c) 2011,2014 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

RGI = RefGeneInfo;
n = getChromosomeNum(RGI,ch);
