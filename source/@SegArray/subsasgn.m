%SUBASGN overloaded method for SegArray class

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

function S = subsasgn(S,I,V)

%   S = SUBSASGN(S,I,V)

switch I.type
  case '()'
    [M,N] = size(S.bpts);
    % cast values to the type of the SegArray
    try
        if isa(V,'SegArray')
            V.vals = cast(V.vals,class(S.vals));
        else
            V = cast(V,class(S.vals));
        end
    catch me
        throwAsCaller(me);
    end
    % handle different numbers of subscripts
    switch length(I.subs)
      case 0
        %% impossible indexing
        error('This error shouldn''t happen! (SUBSASGN w/o indices)');
      case 1
        %% indexing via linear index
        idx = I.subs{1};
        % empty index
        if isempty(idx)
            if max(size(V)) <= 1
                % empty index benign with empty or scalar values
                return;
            else
                % empty index w/bigger-than-scalar value causes error
                try
                    S.bpts(idx) = V;
                catch me
                    throwAsCaller(me);
                end
                errNoError(S);
            end
        end
        if ischar(idx) && strcmp(idx,':')
            idx = (1:numel(S.bpts))'; %! TODO optimize for memory usage
        end
        
        % number of values must match number of indexed locations
        if ~isscalar(V) && nnz(idx) ~= numel(V)
            throwAsCaller(MException('MATLAB:index_assign_element_count_mismatch',...
                                [' In an assignment  A(I) = B, the number of elements in B and\n',...
                                'I must be the same.']));
%!            try % let matlab generate error
%!                S.bpts(1,1) = 1:2;
%!            catch me
%!                throwAsCaller(me);
%!            end
%!            errNoError(S);
        end
        
        if isa(idx,'SegArray') && size(idx,1) ~= 1
            % special-case logical SegArray linear index
            if islogical(idx.vals)
                % empty index leaves array unchanged
                if nnz(idx) == 0
                    return
                end
                % error if value is empty and neither index is
                if isempty(V)
                    throwAsCaller(MException('MATLAB:SegArray:unsupportedFeature', ...
                           'Deleting array elements by assigning empty values is unsupported.'));
                end
                % ensure count of elements == indexed spaces
                if isscalar(V) || sum(idx(:))== numel(V)
                    if ~isscalar(V) && ~isa(V,'SegArray')
                        V = SegArray(V);
                    end
                    
                    S = triple_seg_asgn(S,idx,V);
                    S = anneal(S);
                    return;
                else % generate a "numel must be the same" error
                    try
                        S.bpts(I.bpts) = full(V);
                    catch me
                        throwAsCaller(me);
                    end
                    errNoError(S);
                end
            else
                % to do: break into slices to conserve memory?
                segwarn(idx,'converting linear SegArray index to full for SUBSASGN');
                idx = full(idx);
            end
        end
        if isa(V, 'SegArray')
            % value to assign is a SegArray, but index is not
            segwarn(V,'converting assigned values from SegArray to full');
            V = full(V);
        end
        % special case indexes
        if ischar(I.subs{1})
            if strcmp(I.subs{1},':')
                % ':' case - values must match size or be empty
                if numel(S.bpts) == numel(V)
                    S = reshape(SegArray(V),size(S));
                elseif numel(S.bpts) > 0
                    if isscalar(V)
                        % repeat scalar value V in shape of S
                        S = SegArray.constant(V,M,N);
                    else
                        try % cause element-count-mismatch exception
                            S.bpts(:)=V;
                        catch me
                            throwAsCaller(me);
                        end
                        errNoError(S);
                    end
                end
                return
            else
                idx = cast(idx,'uint8');
            end
        elseif islogical(idx)
            % indexing is logical e.g. from M(M>thresh)
            idx = find(idx);
            if isempty(idx)
                return
            end
        end
        if ~isnumeric(idx)
            if any(strcmp('subsindex',methods(idx))) 
                 idx = subsindex(idx);
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
            throwAsCaller(badSubsException(inputname(1),idx));
        end
        % range test index
        if (0)%(max(idx(:)) > numel(S.bpts)) || (min(idx(:)) < 1)
            % let matlab generate error
            try
                S.bpts(idx);
            catch me
                throwAsCaller(me);
            end
            errNoError(S);
        end

        % indexing is via a numeric full vector
        Jsubs = idx(:);
        biggest = max(Jsubs);
        if biggest > numel(S.bpts)
            % The expansion behavior with one subscript is poorly defined.
            % Vectors seem to be expanded, arrays are not unless they are
            % logical. We convert arrays and row vectors to full and mirror
            % any error.
            if ~isvector(S.bpts) || size(S.bpts,1) == 1
                S = full(S);
                try
                    S(idx) = V;
                catch me
                    throwAsCaller(me);
                end
                % if we survive, convert back to SegArray and exit
                S = SegArray(S);
                return;
            else
                % expand column vector to accommodate largest index
                S = segexpand(S,biggest,1);
            end
        else
            S = indexwrite(S,Jsubs,V(:));
        end
        
      case 2
        %% indexing via two vector subscripts

        idx1 = I.subs{1};
        idx2 = I.subs{2};

        %% process first index (may be logical SegArray, else must be converted to numeric)
        if ischar(idx1)
            if strcmp(idx1,':')
                idx1 = SegArray.constant(true,size(S,1),1);
            else
                idx1 = cast(idx1,'uint8');
            end
        end
        % find maximum index along first dimension
        if islogical(idx1)
            if ~isa(idx1,'SegArray')
                idx1 = find(idx1);
            end
        elseif ~isnumeric(idx1)
            % if not numeric or logical, try to get an index out of it
            if any(strcmp('subsindex',methods(idx1))) 
                idx1 = 1 + subsindex(idx1);
            else
                throwAsCaller(MException('MATLAB:UndefinedFunction',...
                                         'Function ''subsindex'' is not defined for values of class ''%s''.',...
                                         class(idx1) ));
            end        
        end
        % reshape 1st subscript to column vector
        idx1 = subsref(idx1,substruct('()',{':'}));
        % validate values in 1st subscript
        if ~isa(idx1,'SegArray')
            try
                validateattributes(idx1,{'numeric'},{'positive','integer'});
            catch me
                throwAsCaller(badSubsException(inputname(1),idx1));
            end
        end

        %% process second subscript (must be converted to numeric)
        if isa(idx2,'SegArray')
            segwarn(idx2,'converting second SegArray subscript to full for SUBSASGN');
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
        elseif ~isnumeric(idx2)
            if any(strcmp('subsindex',methods(idx2))) 
                idx2 = 1 + subsindex(idx2);
            else
                throwAsCaller(MException('MATLAB:UndefinedFunction',...
                                         'Function ''subsindex'' is not defined for values of class ''%s''.',...
                                         class(idx2) ));
            end        
        end
        % reshape 2nd subscript to column vector
        idx2 = subsref(idx2,substruct('()',{':'}));

        % validate values in 2nd subscript
        try
            validateattributes(idx2,{'numeric'},{'positive','integer'});
        catch me
            throwAsCaller(badSubsException(inputname(1),idx2));
        end

        % if either index is empty, return unchanged array
        if isempty(idx1) || isempty(idx2)
            return
        end
        % error if value is empty and neither index is
        if isempty(V)
            throwAsCaller(MException('MATLAB:SegArray:unsupportedFeature', ...
                    'Deleting array elements by assigning empty values is unsupported.'));
        end
        
        % test that number of values == number of indexed positions
        if ~isscalar(V) && nnz(idx1) * nnz(idx2) ~= numel(V)
            throwAsCaller(MException('MATLAB:subsassigndimmismatch',...
                                     'Subscripted assignment dimension mismatch.'));
        end
        
        % max1 = maximum index along first dimension
        if islogical(idx1)
            if isa(idx1,'SegArray')
                vals = getvals(idx1);
                lnzvalx = find(vals,1,'last');
                if isempty(lnzvalx)
                    max1 = [];
                else
                    if lnzvalx == length(vals)
                        max1 = numel(idx1);
                    else
                        [bpts,~] = find(getbpts(idx1));
                        max1 = bpts(lnzvalx+1)-1;
                    end
                end
            else
                max1 = find(any(idx1,2),1,'last');
            end
        else
            max1 = max(idx1);
        end
        max2 = max(idx2);
        
        % if either subscript exceeds array size, expand array to accommodate
        if ~isempty(max1) && (max1 > M || max2 > N)
            S = segexpand(S,max1,max2);
            [M,N] = size(S);
        end
        
        % special-case logical SegArray as first index
        if isa(idx1,'SegArray')
            if islogical(idx1.vals)
                % if "empty" logical index, return array unchanged
                if ~any(idx1.vals)
                    return
                end
                % error if value is empty and neither index is
                if isempty(V)
                    throwAsCaller(MException('MATLAB:SegArray:unsupportedFeature', ...
                            'Deleting array elements by assigning empty values is unsupported.'));
                end
                % truncate falses at end
                if length(idx1) > M
                    idx1 = subsref(idx1,substruct('()',{1:M,':'}));
                end
                % convert full non-scalar values to SegArray
                if ~isscalar(V) && ~isa(V,'SegArray')
                    V = SegArray(V);
                end
                % TODO below could be more memory efficient, e.g. eliminate idx
                if ischar(I.subs{2}) && strcmp(I.subs{2},':')
                    % no column selection
                    idx = repmat(idx1,1,N);
                    S = anneal(triple_seg_asgn(S,idx,V));
                else
                    % select indexed columns from target SegArray 
                    Scols = subsref(S,substruct('()',{':',idx2}));
                    % assign values to selected columns with repeated SegArray index
                    idx = repmat(idx1,1,size(Scols,2));
                    Scols = anneal(triple_seg_asgn(Scols,idx,V));
                    % now put modified columns back into target SegArray 
                    S = sasgnCols(S,idx2,Scols);
                end
                return;
            else
                % convert numeric SegArray index to full
                segwarn(idx1,'converting first SegArray subscript to full for SUBSASGN');
                idx1 = full(idx1);
           end
        end
        
        % optimized column-level write for first subscript == ':'
        if ischar(I.subs{1}) && strcmp(I.subs{1},':')
            S = sasgnCols(S,idx2,V);
            return
        end

        % create linear index from subscripts
        m = numel(idx1);
        n = numel(idx2);
        if m * n > 0
            idx1 = idx1(:);
            idx2 = reshape(idx2, 1, n);

            % construct linear indices from the two-axes indices
            colsper = max(1,floor(S.MAX_INDICES/m));
            if colsper > N
                % few enough linear indices to do in one pass 
                Jsubs = repmat((idx2-1)*M,m,1) + repmat(idx1,1,n);
                val = full(subsref(V,substruct('()',{':'})));
                S = indexwrite(S,Jsubs(:),val); %!
            else
                % large number of linear indices: do indexing operation in column slices 
                for k0 = 1:colsper:n
                    k1 = min(n,k0+colsper-1);
                    Jsubs = repmat((idx2(k0:k1)-1)*M,m,1) + repmat(idx1,1,k1-k0+1);
                    if isscalar(V)
                        S = indexwrite( S, Jsubs(:), full(V) );
                    else
                        val = full(subsref(V,substruct('()',{':',k0:k1})));
                        S = indexwrite(S, Jsubs(:),val);
                    end
                end
            end
        end
      otherwise
        % TODO index higher than two dimensions???
        throwAsCaller(MException('MATLAB:SEGARRAY:dimTooHighMan', ...
                  'Indexing higher than two dimensions currently unsupported'));
    end
  %% cell and struct access (useless)
  case '{}'
    try % generate error for dereferencing non-cell array
        S.bpts = S{I.subs{:}};
    catch me
        throwAsCaller(me)
    end
  case '.'
    % this destroys the SegArray and makes a struct (what matlab does)
    S = setfield(S,I.subs{1},V);
end

%% Internal helper function for SegArray implementation of SUBSASGN
function S = indexwrite(S,Jsubs,vals)
% S     SEGARRAY being index-written
% Jsubs a column vector of linear indices
% vals  a column vector of values or scalar value
% returns an updated SEGARRAY

% TODO: (1) support SegArray vals argument; 
%       (2) compact duplicate values, especially for scalar vals

% co-sort breakpoints and assignment indices 
[Jsubs,jidx] = sort(Jsubs);            % sort out-of-order indices
if ~isscalar(vals)
    values = [S.vals;vals(jidx)];      % combine values
else
    values = [S.vals;vals];            % put scalar value at end
end
Jarray = find(S.bpts);                 % Jarray linearly indexes nz sparse elements
[sorted,order] = sort([Jarray(:);Jsubs]);
bmask = ( order <= length(Jarray) );   % mark Jarray (breakpoint) indices

% expand a space after each block of contiguous indices for "rebreaking"
igaps = [diff(Jsubs)>1;...             % mark "index gaps" in Jsubs space
         Jsubs(end)<numel(S.bpts)];
xgaps = false(size(bmask));            % expand gaps to co-sort space 
xgaps(~bmask) = igaps;

% create indexed mapping to gap-expanded space
gxlen = length(xgaps) + sum(igaps);
steps = circshift(xgaps,1);
steps(1) = 0;
gxmap = (1:length(bmask))' + cumsum(steps);
% create value index
xbmask = zeros(gxlen,1);
xbmask(gxmap) = bmask;
vidx = cumsum(xbmask);                 % extend breakpoints
if ~isscalar(vals)
    vidx(gxmap(~bmask)) = order(~bmask);
else
    vidx(gxmap(~bmask)) = length(values);
end        
% create position index and fill gaps w/rebreaks
bidx = zeros(gxlen,1);
bidx(gxmap) = sorted;
bidx(~bidx) = sorted(xgaps)+1;
% eliminate all but last of identical adjacent indices
undup = [logical(diff(bidx));true];
bidx = bidx(undup);

% write result values and location indices
S.vals = values(vidx(undup));
[m,n] = size(S.bpts);
[i,j] = ind2sub([m n],bidx);
S.bpts = sparse(i,j,1:numel(S.vals),m,n);


%% optimal columnwise assignment
% S is the SegArray to be updated
% idx2 is a numeric or logical column index
% V is an array of values 
function S = sasgnCols(S,idx2,V)

% force SegArray values, full column index
if ~isa(V,'SegArray')
    V = SegArray(V);
end
if isa(idx2,'SegArray')
    idx2 = full(idx2);
end
    
[Vi,Vj] = find(V.bpts);
allvalues = [S.vals;V.vals];
Sn = length(S.vals);
Vn = length(V.vals);
[P,Q] = size(V.bpts);
% assign columns of a reindexed sparse array of values 
S.bpts(:,idx2) = sparse(Vi,Vj,(Sn+1):(Sn+Vn),P,Q);
% extract merged breakpoints and reordered value index
[I,J,vx] = find(S.bpts);
[M,N] = size(S.bpts);
S.bpts = sparse(I,J,1:length(vx),M,N);
S.vals = allvalues(vx);

%% optimized logical segarray index assignment
function S = triple_seg_asgn(S,I,V)
% S     SEGARRAY being index-written
% I     A logical SegArray index
% V     A SegArray of values or scalar value
% returns an updated SEGARRAY


% scalar value is converted into a SegArray
if isscalar(V)
    V = SegArray.constant(V,nnz(I),1); %! preserves logical
%!  V = V + SegArray.zeros(nnz(I),1);
end

% expand V breakpoints with I
ibrks = find(I.bpts);
ivals = I.vals;
% remove any zero-valued breakpoints at end of I
while ibrks(end) > numel(S) && ~ivals(end)
    % remove 
    ibrks(end) = [];
    ivals(end) = [];
end
bases = ibrks(ivals);
seglens = diff([ibrks;numel(I.bpts)+1]);
iprog = cumsum([1;seglens(ivals)]);
iprog(end)=[];
bases = bases-iprog;
vbrks = find(V.bpts);
[~,order] = sort([iprog;vbrks(:)]);
vmask = order > length(iprog);
iiprog = cumsum(~vmask);
vbrks = vbrks(:)+bases(iiprog(vmask));

% cosort S, I, and expanded V breakpoints
sbrks = find(S.bpts);
[sorted,order] = sort([sbrks(:);vbrks(:);ibrks(:)]);
uniq = [diff(sorted)~=0;true];

% create value indices for each breakpoint
%! TODO could be more memory efficient here
bmask = ( order <= length(sbrks) );   % mark sbrks (breakpoint) indices
istart = length(sbrks) + length(vbrks) + 1;
vmask = (order < istart) & ~bmask;
imask = ( order >= istart);

bidx = cumsum(bmask);
iidx = cumsum(imask);
vidx = cumsum(vmask);

iidx = iidx(uniq);
bidx = bidx(uniq);
vidx = vidx(uniq);

% combine old and new values according to logical index
newflag = ivals(iidx);
values = repmat(V.vals(1),size(newflag)); %! do not force "double" class
values(newflag) = V.vals(vidx(newflag));
values(~newflag) = S.vals(bidx(~newflag));
S.vals = values;

% create new breakpoints
sorted = sorted(uniq);
[m,n] = size(S.bpts);
i = mod(sorted-1,m) + 1;
j = floor((sorted-1)/m) + 1;
S.bpts = sparse(i,j,1:numel(S.vals),m,n);

%% subfunction for expanding a SegArray
% (this is for Matlab's pathological behavior of expanding array size to
% accomodate out-of-range indices)
function S = segexpand(S,max1,max2)
[M,N] = size(S.bpts);

% value to fill margins with
if islogical(S.vals)
    fill = false;
else
    fill = 0;
end
if max1 > M
    % vertcat bottom margin 
    S = [S;SegArray.constant(fill,max1-M,N)];
end
if max2 > N
    % add column header breaks to make right margin
    S.vals = [S.vals; repmat(fill,max2-N,1)];
    S.bpts(1,N+1:max2) = fill;
end
        

%% helper function for bad subscript exception 
function me = badSubsException(name,~)
    access_str = sprintf('index to %s must be a positive integer or logical',name); 
    me = MException('MATLAB:badsubscript',access_str);



