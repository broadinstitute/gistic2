function scores = score_ziggs_by_table(Q,hd,xamp,ylen)
%SCORE_ZIGGS_BY_TABLE look up segment scores from anplitude x length table
%
%   SCORES = SCORE_ZIGGS_BY_TABLE(Q,HD,XAMP,YLEN)
%
% Returns a vector of SCORES corresponding to each row of the ziggurat 
% segments in array Q. HD is a score lookup table, and XAMP and YLEN are
% respectively amplitude and length values for the rows and columns of HD.
%

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

  scores = zeros(size(Q,1),1);
  
  for i=1:size(Q,1)
    % convert amplitude to index
    xidx = find(xamp < Q(i,4),1,'last');
    if isempty(xidx)
      xidx = 1;
    end
    % convert length to index
    yidx = find(ylen < Q(i,8),1,'last');
    if isempty(yidx)
      yidx = 1;
    end
    % look up score
    scores(i) = hd(xidx,yidx);
  end
  
   
