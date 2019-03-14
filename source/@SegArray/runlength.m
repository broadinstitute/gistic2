function rl=runlength(x,segments)
%RUNLENGTH returns the equal valued segments of data vector x.
%
%   RL = RUNLENGTH(X,SEGMENTS) for a vector X returns N X 3 matrix
%   RL, where N is the number of constant-valued segments of vector
%   X. The three columns of the matrix are (1) start row, (2)
%   inclusive end row, (3) value. 
%
%   Optional input SEGMENTS is a M X 2 matrix with each row of
%   SEGMENTS defining the beginning and end of a segment to analyze
%   for further segmentation. 
%
%   If X is a matrix with multiple columns, RUNLENGTH returns a
%   cell array with the runlength result for each column in the
%   corresponding cell.
%
%   Example:
%           x = [ones(1,10) zeros(1,10) -1.*ones(1,10)]
%           rl = runlength(x)
%

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

if isempty(x)
    rl = [];
else
    % make sure x is a SegArray
    if ~isa(x,'SegArray')
        x = SegArray(x);
    end
    % default segments argument is full range single segment
    if ~exist('segments','var')
        if size(x,1)==1   
            segments = [1 size(x,2)];
        else
            segments = [1 size(x,1)];
        end
    else
        if size(segments,2)==1
            segments = runlength(segments);
        end
    end
  
    if min(size(x))>1
        % behavior for array is to create a cell per column
        rl=cell(1,size(x,2));
        for i = 1:size(x,2)
            rl{i} = dorl(subsref(x,substruct('()',{':',i})),segments);
        end
    else
        rl = dorl(x(:),segments);
    end
end

%% DORL Inner runlength function.
% s is a SegArray column vector
% supsegs is a matrix with columns [start,stop] 
function rl = dorl(s,supsegs)

% join adjacent segments with identical values
s = anneal(s,supsegs(:,1));

rl = [];
bptx = find(s.bpts);

% loop over super segments
for j=1:size(supsegs,1)
    % mark breakpoints exclusively after super segment start
    mask = bptx > supsegs(j,1);
    sval = s.vals( max(1,find(mask,1)-1) );

    if isempty(sval)
        sval = s.vals(sum(~mask));
%!      sval = s.vals(1);
    end
    mask = mask & (bptx <= supsegs(j,2));
    cur_bptx = bptx(mask);
    cur_rl = [[supsegs(j,1);cur_bptx], ...
              [cur_bptx-1;supsegs(j,2)], ...
              [sval;s.vals(mask)]];
    rl = [rl; cur_rl];
end
