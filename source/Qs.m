%Qs - GISTIC structure used to hold ziggurat segments (CNA events)
%
% Qs is a structure used to hold the events resolved by ziggurat
% deconstruction of copy number segments. Qs is usually stored in a
% field in the 'D' copy number data structure after ziggurat
% deconstruction, although some lower level functions operate
% directly on Qs. Qs typically consists of four fields, each a
% 12-column ("Q") array:
%
%   Qs.amp contains events classified as pure amplifications
%   Qs.del contains events classified as pure deletions
%   Qs.aod contains events classified as "amplifications over deletions"
%   Qs.doa contains events classified as "deletions over amplifications"
%
% These numeric Q array fields store each event as a row, and each
% column has a special meaning:
%
% COL  DESCRIPTION
%   1  segment chromosome
%   2  segment start position (SNP units)
%   3  segment end position (SNP units)
%   4  SCNA score
%   5  sample identifier
%   6  starting CN level
%   7  ending CN level
%   8  event length as fraction of chromosome arm
%   9  ziggurat deconstruction score
%  10  arm level for event
%  11  -- UNUSED --
%  12  amplitude

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.
