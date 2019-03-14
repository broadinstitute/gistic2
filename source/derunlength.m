function x=derunlength(rl,poslist)
% DERUNLENGTH ('fromRunlength') creates a data array from a runlength array.
%
%   X = DERUNLENGTH(RL,POSLIST) returns a segmented data array
%   created from the three-column array RL listing the data
%   segments: start row in the first column, end row in column two,
%   and value of segment in column three. The output is a column
%   vector with the defined segments.
%
%   If RL is a cell array, then the output is a two-dimensional
%   array and each cell in RL contains the runlength array for a 
%   corresponding output column.
%
%   If the optional POSLIST argument is given, it provides a vector
%   of output positions, and the input segment values are
%   interpreted as indices into these output positions.
%
%   A runlength array is typically created by the RUNLENGTH
%   function.

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


if ~exist('poslist','var')
  poslist=[];
end

%x=zeros(max(rl(:,2)),1);
if iscell(rl)
  if ~isempty(poslist)
    x=nan(poslist(max(rl{1}(:,2))),length(rl));
  else
    x=nan(max(rl{1}(:,2)),length(rl));
  end    
  for i=1:size(x,2)
    if ~isempty(poslist)
      x(:,i)=derunlength(rl{i},poslist);
    else
      x(:,i)=derunlength(rl{i});
    end
  end
else
  if ~isempty(poslist)
    x=nan(poslist(max(rl(:,2))),1);
  else
    x=nan(max(rl(:,2)),1);
  end
  
  last_end=0;
  if ~isempty(poslist)
    for i=1:size(rl,1)
      x(poslist(rl(i,1)):poslist(rl(i,2)))=rl(i,3);
      if rl(i,1)~=last_end+1
        warning('holes in rl');
      end
      last_end=rl(i,2);
    end
  else
    for i=1:size(rl,1)
      x(rl(i,1):rl(i,2))=rl(i,3);
      if rl(i,1)~=last_end+1
        warning('holes in rl');
      end
      last_end=rl(i,2);
    end    
  end
end

