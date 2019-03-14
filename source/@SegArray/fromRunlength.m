function x = fromRunlength(rl,poslist)
% FROMRUNLENGTH ('fromRunlength') creates a SegArray from a runlength array.
%
%   X = FROMRUNLENGTH(RL,POSLIST) returns a SegArray object created
%   from the three-column array RL listing the data segments: start
%   row in the first column, end row in column two, and value of
%   segment in column three. The output is a column vector with the
%   defined segments.
%
%   If RL is a cell array, then the output is a two-dimensional
%   array and each cell in RL contains the runlength array for a 
%   corresponding output column.
%
%   If the optional POSLIST argument is given, it provides a vector
%   of output positions, and the input segment values are
%   interpreted as indices into these output positions.
%
%   A runlength array is typically created by the RUNLENGTH function.
%   
%   This is the SegArray method equivalent of the DERUNLENGTH
%   function, which outputs a full array. 
%

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
    % accumulate the cell contents in a single matrix
    % along with a corresponding column vector
    crl = repmat(0,0,3);
    cols = [];
    for i=1:size(rl,2)
        rli = rl{i};
        crl = [ crl; rli ];
        cols = [cols; repmat(i,size(rli,1),1)];
    end
    x = derun(crl,cols,poslist);
else
    x = derun(rl,1,poslist);
end

function s = derun(rl,cols,poslist)
if size(rl,1) > 1
    % compatible hole check (worth it?)
    if isscalar(cols)
        Eplus = rl(1:end-1,2) + 1;
        B = rl(2:end,1);
    else
        M = max(rl(:,2));
        N = max(cols);
        Eplus = sub2ind([M N],rl(1:end-1,2),cols(1:end-1)) + 1;
        B = sub2ind([M N],rl(2:end,1),cols(2:end));
    end
    if any(Eplus - B)
        warning('holes in rl');
    end
end
% create SegArray
if isempty(poslist)
    s = SegArray.fromSegments(rl(:,1), rl(:,2), cols, rl(:,3));
else
% (NOTE: potential holes anyway with nontrivial poslist - TODO follow up)
    s = SegArray.fromSegments(poslist(rl(:,1)), poslist(rl(:,2)), cols, rl(:,3));
end
