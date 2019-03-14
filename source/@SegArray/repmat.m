function R = repmat(S,varargin)
%REPMAT - SegArray implementation of REPMAT
% The result is always a SEGARRAY.

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


    % some other SegArray argument got us here
    if ~isa(S,'SegArray')
        % convert other arguments and call base repmat
        for k = 1:size(varargin)
            varargin{k} = full(varargin{k});
        end
        R = repmat(S,varargin{:});
        return
    end
    % process repeat size arguments
    if length(varargin) < 1
        %TO DO error requires at least two inputs
    elseif length(varargin) == 1
        % the size is in a single argument vector 
        v = varargin{1};
        if length(v) == 1
            m = v;
            n = v;
        elseif length(v) == 2;
            m = v(1);
            n = v(2);
        else
            %! TODO too many dims err
        end
    elseif length(varargin) == 2
        % the two dimensional size is in two arguments
        m = varargin{1};
        n = varargin{2};
    else
        %! TODO too many dims err
    end

    [M N] = size(S.bpts);
    
    % special case empty repeat size
    if m*n == 0
        R = SegArray.zeros(m*M,n*N);
        R.vals = cast(R.vals,class(S.vals));
        return
    end
    R = S;
    [I J V ] = find(S.bpts);
    % repeat vertically first
    if m > 1
        if all(I==1)
            % special case row SegArray (or all rows identical)
            R.bpts = sparse(I,J,V,m,N);
        else
            % repeat breakpoints vertically
            [I J V] = find(repmat(S.bpts,m,1));
            R.vals = R.vals(V);
            R.bpts = sparse(I,J,(1:length(I))',M*m,N);
        end
    end
    % repeat horizontally
    if n > 1
        R.vals = repmat(R.vals,n,1);
        [I J] = find(repmat(R.bpts,1,n));
        R.bpts = sparse(I,J,(1:length(I))',M*m,N*n);
    end
