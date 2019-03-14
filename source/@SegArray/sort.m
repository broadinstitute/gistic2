function varargout = sort(S,dim,varargin)
%SORT - SegArray implementation of SORT

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


    if ~exist('dim','var') || isempty(dim)
        if size(S.bpts,1) > 1 || isscalar(S.bpts)
            dim = 1;
        else
            dim = 2;
        end
    end
    R = S;
    [M N] = size(S);
    if nargout > 1
        order = zeros(size(S.bpts));
    end
    if dim==1
        [row,column] = find(S.bpts);
        % loop over columns
        for k=1:size(S,2)
            this_col = column==k;
            these_bpts = row(this_col); 
            segsizes = diff([these_bpts;M+1]);
            % if required for output, create order indices
            if nargout > 1
                ordlets = cell(size(this_col));
                for b = 1:size(these_bpts)
                    ordlets{b} = (these_bpts(b):these_bpts(b)+segsizes(b)-1)';
                end
            end
            % sort the column values
            [vals,ord] = sort(S.vals(this_col),dim,varargin{:});
            R.vals(this_col) = vals;
            % calculate reordered breakpoints
            segsizes = segsizes(ord);
            row(this_col) = cumsum([1;segsizes(1:end-1)]);
            % if required, build order index
            if nargout > 1
                order(:,k) = cell2mat(ordlets(ord));
            end
        end
        R.bpts = sparse(row,column,1:length(column),M,N);
    elseif dim == 2
        % "across the grain" sort: compress array vertically
        % to contain only rows with a breakpoint in any column.
        [i,~] = find(S.bpts);
        brkrows = unique(i);
        M = logical(S.bpts);
        M = M(brkrows,:);
        nrows = length(brkrows);
        % create array of values at each "any" breakpoint and sort it
        M = reshape(cumsum(reshape(M,numel(M),1)),nrows,[]);
        R.bpts = repmat(sparse(brkrows,ones(nrows,1),1:nrows,size(S.bpts,1),1),1,N);
        if nargout > 1
            [R.vals order] = sort(S.vals(M),dim,varargin{:});
        else
            if size(M,1)==1
                R.vals = sort(S.vals(M)',dim,varargin{:});
            else
                R.vals = sort(S.vals(M),dim,varargin{:});
            end
        end
        R.vals = R.vals(:);
        R = anneal(R);
            
    else
         % TODO index higher than two dimensions???
         throwAsCaller(MException('MATLAB:SEGARRAY:dimTooHighMan', ...
                     'Indexing higher than two dimensions currently unsupported'));        
    end
    % stuff results into output
    varargout{1} = R;
    if nargout > 1
        varargout{2} = order;
    end
    
