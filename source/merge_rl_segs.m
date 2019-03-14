function rl = merge_rl_segs(rl,rli,val)
%MERGE_RL_SEGS merge two adjacent segments in a RUNLENGTH
%
%   RL = MERGE_RL_SEGS(RL,RLI,VAL)
%
% Operate on the runlength array RL to replace the two segments
% indexed by the pair of indices in RLI with a single segment with
% the new value VAL. Same interface as REMOVE_RL, but runs faster by
% making the assumption that the removed segments are adjacent.
%

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

if size(rl,2)<4
  rl=[rl zeros(size(rl,1),1)];
end

if iscell(rl)
  for i=1:length(rl)
    if ~isempty(rli{i})
      rl{i}=replace_rl(rl{i},rli{i},val{i});
    end
  end
else
  if rl(rli(1),4) == rl(rli(2),4)
    % replace first indexed segment with merge
    rl(rli(1),1) = min( rl(rli(1),1), rl(rli(2),1) ); % start
    rl(rli(1),2) = max( rl(rli(1),2), rl(rli(2),2) ); % end
    rl(rli(1),3) = val;                               % value
    % eliminate second indexed segment
    rl(rli(2),:) = [];
  else
    error('merging segments across chromosomes!');
  end
end
