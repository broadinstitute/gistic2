% SEGARRAY implementation of VERTCAT (concatenation via [ ; ])
% The result is always a SEGARRAY.

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

function R = vertcat(varargin)
    if nargin == 1
        R = varargin{1};
    else
        roffs = zeros(nargin,1); % row offsets
        nbpts = zeros(nargin,1); % breakpoint counts
        N = size(varargin{1},2);
        M = 0;
        %% preprocess arguments
        for k = 1:nargin
            arg = varargin{k};
            % convert to SegArrays
            if ~isa(arg,'SegArray')
                try
                    arg = SegArray(arg);
                catch me
                    throwAsCaller(me);
                end
                varargin{k} = arg;
            end
            [m n] = size(arg.bpts);
            % rig exception if size mismatch
            if n ~= N
                try
                    [arg.bpts;varargin{1}.bpts];
                    errNoError(arg);
                catch me
                    throwAsCaller(me);
                end
            end
            roffs(k) = M;
            M = M + m;
            nbpts(k) = numel(arg.vals);
        end
        %% move data into sorting arrays
        Nvals = sum(nbpts);
        values = cast(zeros(Nvals,1),class(varargin{1}.vals));
        Jarray = zeros(Nvals,1);
        vi = 0;
        for k = 1:nargin
            S = varargin{k};
            values(vi+1:vi+nbpts(k)) = S.vals;
            [I J] = find(S.bpts);
            Jarray(vi+1:vi+nbpts(k)) = sub2ind([M N],I+roffs(k),J);
            vi = vi + nbpts(k);
        end
        %% sort indices and create output array
        [sorted order] = sort(Jarray);
        R = SegArray;
        R.vals = values(order);
        [I J] = ind2sub([M N],sorted);
        R.bpts = sparse(I,J,1:Nvals,M,N);
    end
end

