function R = subsref(S,I)
% SEGARRAY implementation of SUBSREF (indexed read)

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


% Outer function does argument processing and validation, 
% then calls the appropriate subfunction

    [M N] = size(S.bpts);
    
    switch I.type
      case '()'
        switch length(I.subs)
          case 0
            R = S;
          
          case 1
            % "linear indexing" by one index
            idx = I.subs{1};
            if isempty(idx)
                R = SegArray;
                R.bpts = sparse([]);
                R.vals = cast([],class(S.vals));
                return;
            end
            if isa(idx,'SegArray')
                if size(S.bpts,1) == 1
                    % S is a row vector SegArray (i.e. result of compressing 
                    % a matrix SegArray on dimension 1)
                    try
                        R = S;
                        R.bpts = S.bpts(idx.vals);
                        R.vals = S.vals(idx.vals);
                        catch me
                        throwAsCaller(me);
                    end
                    return;
                else
                    % S is a column vector or matrix
                    if islogical(idx.vals)
                        S.bpts = S.bpts(:);
                        idx = anneal(idx);
                        idx.bpts = idx.bpts(:);
                        % range test index
                        if length(idx.bpts) < length(S.bpts)
                            % pad index shortfall with zeros
                            idx = [idx;SegArray.false(length(S.bpts)-length(idx.bpts),1)];
                        elseif length(idx.bpts) > length(S.bpts)
                            % make sure no nonzeros beyond end of matrix
                            L = length(S.bpts);
                            last1 = find(idx.vals,1,'last');
                            fidx = find(idx.bpts);
                            if isempty(fidx)
                                R = [];
                                return;
                            else
                                % check for NZ breakpoints beyond end of matrix
                                if (fidx(last1) > L)
                                    throwAsCaller(badSubsException);
                                end
                                %make sure last NZ segment ends before end of S
                                if last1 >= length(fidx) || fidx(last1+1) <= L+1
                                    throwAsCaller(badSubsException);
                                end
                            end
                        end
                        R = losexread(S,idx);
                        return
                    else
                        segwarn(idx,'converting linear SegArray index to full for SUBSREF');
                        idx = full(idx);
                    end
                end
            end
            if ischar(idx)
                if strcmp(':',idx)
                    % ':' case - output is reshaped to SegArray column vector
                    R = reshape(S,numel(S),1);
                    return;
                else
                    % convert all other character indices to numeric indices
                    idx = cast(idx,'uint8');
                end
            elseif islogical(idx)
                % logical indexing: convert to numeric index
                idx = find(idx);
                if isempty(idx)
                    R = SegArray;
                    R.bpts = sparse(reshape([],0,1)); % shape matlab returns
                    R.vals = cast([],class(S.vals));
                    return;
                end
            end
            if ~isnumeric(idx)
                if any(strcmp('subsindex',methods(idx))) 
                    idx = subsindex(idx) + 1;
                else
                    throwAsCaller(MException('MATLAB:UndefinedFunction',...
                                             'Function ''subsindex'' is not defined for values of class ''%s''.',...
                                             class(idx) ));
                end        
            end
            % validate numeric index
            try
                validateattributes(idx,{'numeric'},{'positive','integer'});
            catch me
                throwAsCaller(badSubsException);
            end
            % range test index
            if (max(idx(:)) > numel(S.bpts)) ||  (min(idx(:)) < 1)
                try  % let matlab generate error
                    S.bpts(idx);
                catch me
                    throwAsCaller(me);                                
                end
                errNoError(S);
            end
            % nuanced shaping behavior of vector indexing vector:
            % indexed vector (row vs column) wins.
            if isvector(S.bpts) && isvector(idx)
                if 1==size((S.bpts),1) &&  1 < size(idx,1)
                    idx = reshape(idx,1,numel(idx));
                elseif N == 1 && 1 < size(idx,2)
                    idx = idx(:);
                end
            end
            % call inner function (finally!)
            R = indexread(S,idx);

          case 2
            %% indexing via two subscripts
            idx1 = I.subs{1}(:);
            idx2 = I.subs{2}(:);

            %% process first subscript
            if isa(idx1,'SegArray')
                if islogical(idx1)
                    % optimized subfunction for logical SegArray
                    idx1 = anneal(idx1);
                    % -> test index out-of-range
                    last1 = find(idx1.vals,1,'last');
                    if last1 == numel(idx1.vals);
                        last1 = numel(idx1);
                    else
                        bpi = find(idx1.bpts);
                        last1 = bpi(last1+1);
                    end
                    if last1 > M
                        % if index is too large, generate matlab error
                        try
                            S.bpts(last1,:);
                        catch me
                            throwAsCaller(me);
                        end
                        errNoError(S);
                    end
                    if size(idx1.bpts,1) < M
                        % if index is smaller than size, pad to correct size w/falses
                        idx1 = [idx1;SegArray.false(M-size(idx1,1),1)];
                    end
                    % call inner function
                    if isequal(idx2,':')
                        R = losexread(S,idx1);
                    else
                        R = losexread(swizzleCols(S,idx2),idx1);
                    end
                    return;
                else
                    % first subscript is non-logical SegArray
                    segwarn(idx1,'converting first index from SegArray to full for SUBSREF');
                    idx1 = full(idx1);
                end
            end
            if ischar(idx1)
                if strcmp(':',idx1)
                    % use optimized pure column reordering
                    R = swizzleCols(S,idx2);
                    return;
                else
                    idx1 = cast(idx1,'uint8');
                end
            elseif islogical(idx1)
                idx1 = find(idx1);
           %else TODO non-numeric subsindex
            end

            % process second subscript index
            if isa(idx2,'SegArray')
                segwarn(idx2,'converting second index from SegArray to full for SUBSREF');
                idx2 = full(idx2);
            end
            if ischar(idx2)
                if strcmp(':',idx2)
                    idx2 = 1:N;
                else
                    idx2 = cast(idx2,'uint8');
                end
            elseif islogical(idx2)
                idx2 = find(idx2);
           %else TODO non-numeric subsindex
            end

            m = numel(idx1);
            n = numel(idx2);
            % 2D result with a zero dimension 
            if m*n == 0
                R = SegArray;
                R.bpts = sparse([],[],0,m,n);
                R.vals = cast([],class(S.vals));
                return;
            end

            idx1 = idx1(:);
            idx2 = reshape(idx2, 1, n);

            % range test indices and make sure they are positive integers
            if max(idx1) > M         || max(idx2) > N ||... 
               min(idx1) < 1         || min(idx2) < 1 ||...
               any(rem(idx1,1) ~= 0) || any(rem(idx2,1) ~=0)
                try
                    S.bpts(idx1,idx2); % let matlab generate error
                catch me
                    throwAsCaller(me);
                end
                errNoError(S);
            end

            % construct linear index from the two-axis indices
            colsper = max(1,floor(S.MAX_INDICES/m));
            if colsper > N
                % few enough linear indices to do in one pass 
                Jsubs = repmat((idx2-1)*M,m,1) + repmat(idx1,1,n);
                R = indexread(S,Jsubs);
            else
                % large number of linear indices: do indexing operation in column slices 
                R = SegArray.zeros(m,0);
                for k0 = 1:colsper:n
                    k1 = min(n,k0+colsper-1);
                    Jsubs = repmat((idx2(k0:k1)-1)*M,m,1) + repmat(idx1,1,k1-k0+1);
                    R = [ R, indexread(S,Jsubs) ];
                end
            end

          otherwise
             % TODO index higher than two dimensions???
             throwAsCaller(MException('MATLAB:SEGARRAY:dimTooHighMan', ...
                         'Indexing higher than two dimensions currently unsupported'));
        end

      % unsupported indexing types
      case '{}'
        throwAsCaller(MException('MATLAB:cellRefFromNonCell', ...
                                 'Cell contents reference from a non-cell array object.')); 
      case '.'
        throwAsCaller(MException('MATLAB:nonStrucReference', ...
                                 'Attempt to reference field of non-structure array.'));
      otherwise
        throwAsCaller(MException('MATLAB:subsTypeMustBeSquigglyOrSmooth', ...
                         ['The "type" field for the subscript argument to SUBSREF and SUBSASGN\n', ...
                          'must be a character array  of "." or "{}" or "()".'] ));
    end
end

%% INDEXREAD: inner subscripted read operation

function R = indexread(S,Jsubs)
% S     SEGARRAY being index-read
% Jsubs a shaped array of linear indices
% R     the resulting SEGARRAY of indexed values, or a full array
%    if Jsubs is a scalar or empty
%       The shape of R is set to the shape of Jsubs

    % note shape of indices and flatten them to column vectors
    [m n] = size(Jsubs);
    Jarray = find(S.bpts);                    % Jarray linearly indexes nz elements
    if isscalar(Jsubs)
        % index scalar
        idx = find(Jarray <= Jsubs,1,'last'); % find the index before
        R = S.vals(idx);
    else
        % index SegArray
        Jsubs = Jsubs(:);
        jsorted = issorted(Jsubs);
        if (~jsorted)
            [Jsubs,jidx] = sort(Jsubs);       % sort out-of-order indices
        end
        bmask = (SegArray.cosort(Jarray,Jsubs) == 1);  %!
%!      [~,order] = sort([Jarray(:);Jsubs]);
%!      bmask = ( order <= length(Jarray) );  % mark Jarray ("from") indices
        
        bidx = cumsum(bmask);
        % create vector of linear breakpoint indices
        if jsorted
            fuvx = bidx(~bmask);
        else
            fuvx = zeros(m*n,1);
            fuvx(jidx) = bidx(~bmask);        % unsort full array of value indices
        end
        bpmask = [true;logical(diff(fuvx))];  % logical mask of index breakpoints
        bpmask(1:m:(m*n)) = true;             % add column breakpoints
        % put indices and values into SegArray
        R = SegArray();
        N = sum(bpmask);
        [i j] = ind2sub([m n],find(bpmask));  %! reshape breakpoints
        R.bpts = sparse(i,j,1:N,m,n);
        R.vals = S.vals(fuvx(bpmask));
    end
end


%% LOSEXREAD: sparse logical indexed read operation

function R = losexread(S,losex)
% S     SEGARRAY being indexed ('data')
% losex logical SegArray column vector index 
% R     the resulting SEGARRAY of indexed values, or a full array
%
% The shape of R is set to the shape of losex.

    % note shape of indices and flatten them to column vectors
    m = sum(losex);            %! calling SegArray method may be slow
    n = size(S.bpts,2);
    % dispense with the empty result case immediately
    if m * n == 0
        R = SegArray;
        R.bpts = reshape([],m,n);
        R.vals = S.vals([]);
        return;
    end

    % Jarray linearly indexes data breakpoints
    Jarray = find(S.bpts);
    % Larray linearly indexes logical index breakpoints
    Larray = bsxfun(@plus,(0:n-1)*size(S.bpts,1),find(losex.bpts));

    % co-sort locations of data and index breakpoints
%! test cosort method ...
    [group,sorted] = SegArray.cosort(Jarray,Larray); %!
    bmask = (group == 1);                   %!
%!  [sorted,order] = sort([Jarray(:);Larray(:)]);
%!  bmask = (order <= numel(Jarray));   % mark data bpt locations

    bidx = cumsum(bmask);               % numerical index for S.vals
    lidx = cumsum(~bmask);              % numerical index for losex.vals

    mi = numel(losex.vals);
    % get values for logical SegArray subscript
    lmask = losex.vals(mod(lidx-1,mi)+1);
    % remove duplicated position indices
    lmask = lmask & [logical(diff(sorted));true];
    bidx = bidx(lmask);                 % compress value index
    sorted = sorted(lmask);
    % remove all but first of duplicated value indices
    lmask = [true;logical(diff(bidx))];
    bidx = bidx(lmask);
    sorted = sorted(lmask);
    % put values and indices into SegArray
    R = SegArray();
    R.vals = S.vals(bidx);
    N = length(bidx);
    [i j] = ind2sub(size(S.bpts),sorted);

    % map/compress the row index for new breakpoints
    map = cumsum(full(losex)); %!!! note this is slow
    i = map(i);
    R.bpts = sparse(i,j,1:N,m,size(S.bpts,2));
end

%% SWIZZLECOLS: reorder columns only
% Optimize "column swizzle" where all rows are selected and
% no breakpoints need to be synthesized. The algorithm is to
% reorder the sparse bpts array (rethrow any indexing errors),
% and use the new order of nonzeros to reorder the values.
% (LOGBTS would need to synthesize indexing breakpoints if bpts became logical)
function R = swizzleCols(S,idx2)
    if isa(idx2,'SegArray')
        idx2 = full(idx2);
    end
    try
        swizzle = S.bpts(:,idx2);
    catch me
        throwAsCaller(me);
    end
    [i j P] = find(swizzle);
    R = SegArray();
    R.vals = S.vals(P);
    R.bpts = sparse(i,j,1:size(P),size(swizzle,1),size(swizzle,2));
end

function me = badSubsException
    me = MException('MATLAB:badsubscript',...
         'Subscript indices must either be real positive integers or logicals.');
end
