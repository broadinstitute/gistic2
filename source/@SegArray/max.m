function [C,I] = max(A,B,dim)
% SegArray implementation of MAX
%
%    [C,I] = MAX(A,B,dim)
%
%    If A and B are given and at least one is a non-scalar
%    SegArray, a SegArray is returned.
%
%    If just A is given and dim=1, the results are full arrays.  
%    If dim=2, the results are SegArrays.
%

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


% NOTE: maintained in parallel with MIN
    
    % process dimension argument
    if exist('dim','var')
        if ~isempty(B)
            try % generate appropriate error
                max(1,1,1);
            catch me
                throwAsCaller(me);
            end
            errNoError(SegArray);
        else
            try 
                validateattributes(dim,{'numeric'},{'scalar','integer','positive'});
                if dim > length(size(A))
                    throw(MException('dont:care','really dont'));
                end
            catch me
                try % generate appropriate error
                    max(pi,[],full(dim))
                catch me
                    throwAsCaller(me)
                end
            end
        end
    else
        if size(A.bpts,1)==1
            dim = 2;
        else
            dim = 1;
        end
    end

    % process second argument (present if comparing 2 arrays)
    if exist('B','var') && ~isempty(B)
        if nargout == 1
            C = binary(A,B,@max); % TODO clean up exception processing
        else
            try % generate appropriate error
                [C,I] = max(pi,pi);
            catch me
                throwAsCaller(me);
            end
            errNoError(SegArray);
        end
    else

        % max of single array
        if isscalar(A)
            C = A.vals;
            I = 1;
        else
            [M,N] = size(A.bpts);
            [K,J] = find(A.bpts);
            switch dim
              case 1 %% maximum across columns ("with the grain")
                % with the grain
                if isempty(A)
                    % empty
                    if M > 0
                      M = 1;
                    end
                    C = zeros(M,N);
                    I = zeros(M,N);
                elseif M == 1
                    % row vector
                    C = A.vals';
                    I = ones(1,N);
                elseif N == 1
                    % column vectors are special
                    [C I] = max(A.vals);
                    I = K(I);
                else
                    C = accumarray([ones(length(J),1) J],A.vals,[1 N],@max,-Inf);
                    if nargout > 1
                        I = -Inf(size(C));
                        for col=1:length(I)
                            I(col) = K(find(J==col & A.vals==C(col), 1));
                        end
                    end
                end
              case 2 %% maximum across rows ("against the grain")
                if isempty(A)
                    % empty
                    if N > 1
                        N = 1;
                    end
                    C = SegArray.zeros(M,N);
                    I = SegArray.zeros(M,N);
                else
                    [C I] = crossgrain2(A,@(x) max(x,[],2));
                end
              otherwise
                % higher dimensions unsupported
                throwAsCaller(MException('MATLAB:SEGARRAY:dimTooHighMan',...
                                         'Dimensions higher than two unsupported'));
            end
        end
    end
