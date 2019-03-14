function [V I J] = unique(S, varargin)
% UNIQUE - SegArray implementation of the unique function
%   [V I J] = UNIQUE(S,...) emulates built-in UNIQUE
%   V and I are returned as full arrays, J is a SegArray

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


% TODO fix UNIQUE to imitate multiple NaN behavior?
    
    % let built in unique function validate the option args
    try
        unique(pi,varargin{:});
    catch me
        throwAsCaller(me);
    end
    if size(S.bpts,1) == 1
        [V I J] = unique(full(S),varargin{:});
    else
        % evaluate options
        first_occurrence = false;
        by_row = false;
        for i = 1:length(varargin)
            if strcmp(varargin{i},'first')
                first_occurrence = true;
            elseif strcmp(varargin{i},'rows')
                by_row = true;
            end
        end
        J = SegArray;
        if by_row
            % operate on reduced array of rows with breakpoints
            [allrows,~] = find(S.bpts);
            brkrows = unique(allrows);
            brkS = subsref(S,substruct('()',{brkrows,':'}));
            [V Iu Ju] = unique(full(brkS),varargin{:});
            if first_occurrence
                I = brkrows(Iu);
            else % last occurrance
                brkrows_ = [brkrows;length(brkrows)+1];
                I = brkrows_(Iu+1)-1;
            end
            J.bpts = brkrows;
        else % by element
            [V Iu Ju] = unique(S.vals,varargin{:});
            K = find(S.bpts);
            if first_occurrence
                I = K(Iu);
            else % last occurrance
                K = [K;numel(S.bpts)+1];
                I = K(Iu+1)-1;
            end
            J.bpts = S.bpts(:);
        end
        J.vals = Ju;
    end
