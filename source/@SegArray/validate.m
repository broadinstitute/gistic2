function validate(S,saname)
% VALIDATE  internal SegArray function for validating SegArray
%     VALIDATE(S,SANAME) does nothing if S is a valid SegArray, or throws
%     an exception if S does not have a valid internal representation.
%     SANAME is an optional string to use in the exception message
%     for referring to the SegArray variable, the default is the variable 
%     name in the caller's workspace.

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


    if ~exist('saname','var') || isempty(saname)
        saname = inputname(1);
    end
    % S.bpts must be sparse
    if ~issparse(S.bpts)
        throwAsCaller(MException('MATLAB:SegArray:internalRepError',...
                                 '%s.bpts should be a sparse array',...
                                 saname));
    end
    % S.vals should not be sparse
    if issparse(S.vals)
        throwAsCaller(MException('MATLAB:SegArray:internalRepError',...
                                 '%s.vals should not be a sparse array',...
                                 saname));
    end
    % number of values must equal number of breakpoints (critical)
    if nnz(S.bpts) ~= numel(S.vals)
        throwAsCaller(MException('MATLAB:SegArray:internalRepError',...
                                 'Number of breakpoints should equal number of values'));
    end
    % non-empty S.vals must be a column vector
    if ~isempty(S.vals) && size(S.vals,2) ~= 1
        throwAsCaller(MException('MATLAB:SegArray:internalRepError',...
                                 '%s.vals should be a column vector',...
                                 saname));
    end
    % validate breakpoint attributes (non-negative integers)
    try
        validateattributes(S.bpts,{'numeric'},{'integer','nonnegative','nonnan'});
    catch me
        throwAsCaller(MException('MATLAB:SegArray:internalRepError',...
                                 '%s.bpts attribute error: %s',...
                                 saname,me.message));
    end
    % validate value class
    if ~isnumeric(S.vals) && ~ischar(S.vals) && ~islogical(S.vals)
        throwAsCaller(MException('MATLAB:SegArray:internalRepError',...
                                 '%s.vals should be numeric, logical, or character',...
                                 saname));
    end
        
        
