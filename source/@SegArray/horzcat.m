% SEGARRAY implementation of HORZCAT (concatenation via [ ])
% The result is always a SEGARRAY.

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

function R = horzcat(varargin)
    if nargin == 1
        R = varargin{1};
    else
        coffs = zeros(nargin,1); % column offsets
        nbpts = zeros(nargin,1); % breakpoint counts
        M = size(varargin{1},1);
        N = 0;
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
            if m ~= M
                try
                    horzcat(arg.bpts,varargin{1}.bpts);
                    errNoError(arg);
                catch me
                    throwAsCaller(me);
                end
            end
            coffs(k) = N;
            N = N + n;
            nbpts(k) = numel(arg.vals);
        end
        %% Create SegArray and directly move each array input to it
        % (no sort required, arrays are already in output order)
        Nvals = sum(nbpts);
        R = SegArray;
        R.vals = zeros(Nvals,1);
        R.bpts = sparse([],[],[],M,N,Nvals);
        vi = 0;
        for k = 1:nargin
            S = varargin{k};
            R.vals(vi+1:vi+nbpts(k)) = S.vals;
            I = find(S.bpts);
            R.bpts( I + M*coffs(k) ) = vi+1:vi+nbpts(k);
            vi = vi + nbpts(k);
        end
    end
end

