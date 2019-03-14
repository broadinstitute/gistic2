classdef (InferiorClasses = {?double,?single,?logical,?char}) SegArray 
%SEGARRAY Create compressed segmented array from normal array.
%   S = SEGARRAY(M) converts a two-dimensional matrix M into a
%   SegArray representation, which efficiently stores matrix data which
%   contains repeats of the same value along the first dimension.
%
%   S = SEGARRAY creates an empty SegArray.
%  
%   SegArrays should behave like ordinary arrays within the scope of
%   their usage by GISTIC.
%
%   See also SegArray.fromSegments, SegArray.fromRunlength, SegArray.zeros
%   for other ways to make SegArrays.

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


%  4 Mar 08 -- jdobson@broad.mit.edu
% $HeadURL: https://svnrepos/CancerGenomeAnalysis/trunk/gistic2core/source/@SegArray/SegArray.m $
% $Id: SegArray.m 82894 2016-04-14 16:29:15Z schum $
% $Rev: 82894 $
% $LastChangedBy: schum $ @broadinstitute.org
% $Date: 2016-04-14 12:29:15 -0400 (Thu, 14 Apr 2016) $

    properties (SetAccess = 'public', GetAccess = 'public')
        vals = [];  % breakpoint values and their class
        bpts = [];  % array shape and breakpoint locations
    end
    properties (Constant = true)
        % MAXINDICES throttles memory usage by subsref and subsasgn.
        % Indexed columns are grouped so that the number of indices is 
        % below MAX_INDICES.
        MAX_INDICES = 1e7;
    end
    methods (Access = 'public')
        %% constructor
        function S = SegArray(B)

            if nargin == 0
                % create empty SegArray
                S.vals = [];
                S.bpts = [];
            elseif nargin == 1
                % create from existing array
                if isa(B,'SegArray')
                    % done already
                    S = anneal(B);
                    return;
                end
                B = full(B);
                [M N] = size(B);
                if isempty(B)
                    S.bpts = sparse(M,N);
                    S.vals = reshape(B,0,0); % keep class in vals
                elseif M == 1
                    S.vals = B(:);
                    S.bpts = sparse(1,1:N,1:N,M,N);
                else
                    % find segment breakpoints, including NaN segments
                    Xchanges = diff(B) ~= 0;
                    Xchanges = [true(1,N);Xchanges];
                    nanm = isnan(B);
                    Xchanges(nanm & circshift(nanm,1)) = false;
                    Xchanges(1,:) = true;            % column breakpoints

                    % create sparse array of breakpoint value indexes
                    [I J] = find(Xchanges);
                    S.bpts = sparse(I,J,1:length(I),M,N);

                    % store values in SegArray
                    Xchanges = B(Xchanges);
                    S.vals = Xchanges(:);
                end
            else
                % too many arguments
                throwAsCaller(MException('MATLAB:SEGARRAY:tooManyArgs',...
                          'Too many arguments for SEGARRAY constructor.'));
            end
        end
    end
    %% public static methods
    methods (Static)
        % factory method for making a SegArray from segment descriptors
        S = fromSegments(B,E,C,V,F,M,N,TOL);
        % factory method for making a SegArray from a run-length array
        x = fromRunlength(rl,poslist);
        % factory method for making an all-zeros MxN SegArray
        S = zeros(M,N);
        % factory method for making an all-false MxN SegArray
        S = false(M,N);
        % factory method for making a constant SegArray
        S = constant(val,M,N);

        % return revision and version numbers
        vertext = version;
    end
    methods (Access = 'public')
        % overloaded subsref/subsasgn/display functions
        
        b = subsref(S,s);

        display(S);

        disp(S);
        
        S = subsasgn(S,s,b);
        
        I = subsindex(A);

        % other essential overloads (%!ses)

        M = full(S); % convert to full array
        
        %% class-specific interface methods

        % get breakpoints in a sparse array
        bpts = getbpts(S);
        % get breakpoint counts (per column)
        nbpts = getbpt_counts(S);
        % get values in a column vector
        vals = getvals(S);
        % get ratio of represented values to breakpoints
        compression_factor = get_compression_factor(S);
        % get number of values (redundant?)
        num_vals = get_numvals(S);
        % join adjacent segments with the same value
        S = anneal(S,segments);
        % get a runlength representation of the SegArray
        rl = runlength(x,segments);
        % sum segment list into numeric SegArray
        S = addSegments(S,B,E,C,V);
        % get segment list from SegArray
        [B,E,C,V] = getSegments(S,F);
        % limit values to specified range
        S = cap_vals(S,caps);
        % validate internal representation of SegArray
        validate(S,varname);
 
        %% Other functions

        R = reshape(S,varargin);
        sumvect = sum(S,dim);
        sumvect = nansum(S,dim);
        nz = nnz(S);
        R = repmat(S,varargin);
        
        C = mtimes(A,B);
        
        C = bsxfun(fun,A,B);
        
        %% binary operations

        % private helper function for simple array-wise or scalar operations 
        C = binary(A,B,op); % TODO make private
        
        C = plus(A,B);
        C = minus(A,B);
        C = ldivide(A,B);
        C = rdivide(A,B);
        C = mod(A,B);
        C = power(A,B);
        C = gt(A,B);
        C = lt(A,B);
        C = eq(A,B);
        C = ge(A,B);
        C = le(A,B);
        C = ne(A,B);
        C = and(A,B);
        C = or(A,B);
        C = xor(A,B);

        %% arraywise unary operations

        B = uminus(A);
        TF = not(A);
        B = log(A);
        B = log2(A);
        B = log10(A);
        B = exp(A);
        B = reallog(A);
        TF = logical(A);
        TF = isnan(A);
        TF = isinf(A);
        B = floor(A);
        B = ceil(A);
        B = round(A);
        B = fix(A);
        B = double(A);
        B = real(A);
        B = imag(A);
        B = abs(A);
        B = sign(A);
        B = sqrt(A);

        %% data queries

        TF = isnumeric(A);
        TF = islogical(A);
        TF = ischar(A);
        TF = isreal(A);

        %% shape queries

        varargout = size(S,varargin);
        varargout = numel(S,varargin);
        L = length(S);
        TF = isscalar(A);
        TF = isvector(A);
        TF = isempty(A);

        %% complex operations
        
        R = crossgrain(S,ccfh); %TODO make private
        [Q R] = crossgrain2(S,ccfh); %TODO make private
        
        Q = mrdivide(A,B); % exists for poorly written built-ins like median
        B = any(A,dim);
        B = all(A,dim);
        varargout = unique(varargin);
        [I,J,V] = find(X,K,str); % TODO more testing, fix NaN values
        varargout = sort(varargin);
        len = end(S,k,n);
        y = mean(x,dim);
        y = nanmean(x,dim);
        Y = diff(X,n,dim);
        varargout = min(varargin);
        varargout = max(varargin);
        S = cat(DIM,varargin);
        S = vertcat(varargin);
        S = horzcat(varargin);
        
        meds = nanmedian(S,dim);
        
        %% pass-through functions

        varargout = histc(varargin);
%        varargout = image(varargin);
%        varargout = imagesc(varargin);

    end
    %% private methods
    methods (Access = 'private')
        S = errNoError(S);  % if rigged exception doesn't work
        segwarn(S,msg); % TODO this is slow, even when messages are off
    end
    %% helper Mex methods (now public)
    methods (Static,Access = 'public')
        [g,v,i] = cosort(varargin);
        values = sumranges(B,E,V,valsize);
    end
end
