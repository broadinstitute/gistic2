function snps=find_snps(C,chrn,st,en,flanking)
% FIND_SNPS returns marker positions for a single genomic range.
%
%   SNPS = FIND_SNPS(D,CHRN,START,END,FLANKING) returns a column
%   vector of marker base positions corresponding to the range
%   [START,END] on chromosome CHRN. 
%
%   If FLANKING=0 (default), then only markers included in the
%   range [START,END] are returned. If FLANKING=1, ranges without
%   markers will append the markers immediately above and below the
%   range. If flanking=2, flanking markers will be added even if
%   markers have been found in the range.
%
%   All arguments except D must be scalars.

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


% ---
% $Id$
% $Date: 2007-12-03 15:57:57 -0500 (Mon, 03 Dec 2007) $
% $LastChangedBy: rameen $
% $Rev$

if ~isfield(C,'chrn')
  C=add_chrn(C);
end
if ~exist('flanking','var')
  flanking=0;
end

if is_mb(C)
  C.pos=C.pos*1e6;
  warning('Dont use MB anymore.');
end

if ischar(chrn)
  tmp=regexp(chrn,'chr([0-9XY]+):([0-9]+)-([0-9]+)','tokens');
  chrn=chromosome2num(tmp{1}{1});
  st=str2num(tmp{1}{2});
  en=str2num(tmp{1}{3});
end

in_chr=find(C.chrn==chrn);
if isempty(in_chr)
% warning('empty chromosome');
  snps=[];
  return
end

cpos = C.pos(in_chr);
snps = find(cpos>=st & cpos<=en);
snps = snps';

if isempty(snps) && flanking
  snps1=find(cpos<st);
  if ~isempty(snps1)
    snps(1)=snps1(end);
  else
%!    warning('no snp before');
    snps(1)=1;
  end
  snps2=find(cpos>en);
  if ~isempty(snps2)
    snps(2)=snps2(1);
  else
%!   warning('no snp after');
    snps(2)=length(in_chr);
  end
  if isempty(snps)
    disp('still empty');
  end
elseif ~isempty(snps) && flanking==2
  snps1=find(cpos<st);
  if ~isempty(snps1)
    snps=[ snps1(end) snps];
  end

  snps2=find(cpos>en);
  if ~isempty(snps2)
    snps=[ snps snps2(1)];
  end
end

if ~isempty(snps)
  snps=in_chr(snps);
end
