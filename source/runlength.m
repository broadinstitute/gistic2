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
%   for further segmentation. If SEGMENTS is a column vector, it
%   is first converted into a runlength array describing its
%   segments. 
%
%   If X is a matrix with multiple columns, RUNLENGTH returns a
%   cell array with the runlength result for each column in the
%   corresponding cell.
%
%   DERUNLENGTH provides the inverse functionality.
%
%   Example:
%           x = [ones(1,10) zeros(1,10) -1.*ones(1,10)]
%           rl = runlength(x)
%

% GISTIC software version 2.0
% Copyright (c) 2011, 2016 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


if isempty(x)
    rl = [];
else
    if ~exist('segments','var')
        if size(x,1)==1   
            segments = [1 size(x,2) 1];
        else
            segments = [1 size(x,1) 1];
        end
    else
        if size(segments,2)==1
            segments = runlength(segments);
        end
    end
  
    if min(size(x))>1
        rl = cell(1,size(x,2));
        for i=1:size(x,2)
        %     disp(i);
            rl{i} = runlength(x(:,i),segments);
        end
    else
        x=as_row(x);
        rl=[];
        for j=1:size(segments,1)  %loop over rows of segments vector  (1)

            y = x(segments(j,1):segments(j,2));
            [u,~,uj] = lunique(y);
            %u(uj) = y for sorted u, so uj gives rank of each element in y

            % combine NaNs to segments
            nanpos = find(isnan(u));
            if length(nanpos)>1
                % fix only uj
                uj(isnan(y)) = max(uj(isnan(y)));
            end
            ch = find(abs(diff([0 uj 0])));
            % ch(i) is 1 if uj(i+1)-uj(i)~=0 (if input x is 1,0,-1 for pos,
            % zero and neg slope, ch gives indices of change in x)

            cur_rl = zeros(length(ch)-1,3);
            cur_rl(:,1) = ch(1:(end-1));  %start of segment
            cur_rl(:,3) = y(cur_rl(:,1)); %values for segment
            cur_rl(:,1) = cur_rl(:,1)+segments(j,1)-1;
            cur_rl(:,2) = ch(2:end)-1+segments(j,1)-1; %end of the segment
            rl = [rl; cur_rl];
        end % for loop over segments
    end
    %diff gives 1st difference so that if dx = diff(x), dx(n) = x(n+1)-x(n)
end
