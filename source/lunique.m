function varargout = lunique(varargin)
%LUNIQUE - unique() wrapper to deal with the version lunacy of the 'legacy' option

% GISTIC software version 2.0
% Copyright (c) 2011-2017 Gad Getz, Rameen Beroukhim, Craig Mermel,
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, Gordon Saksena
% All Rights Reserved.
% (See the accompanying LICENSE file for licensing details.)

% check matlab version to see if legacy option is supported
if str2double(regexprep(version,'\.[0-9]+\.[0-9]+ .+$','')) >= 7.14
    % make sure it's not already an option
    noleg = true;
    for i=1:length(varargin)
        if strcmp(class(varargin{i}),'char') && strcmp('legacy',varargin{i})
            noleg = false;
            break;
        end
    end
    if noleg
        % add 'legacy' option
        varargin = [varargin,'legacy'];
    end
end
[varargout{1:nargout}] = unique(varargin{:});