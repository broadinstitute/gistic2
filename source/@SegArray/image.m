function varargout = image(varargin)
% SegArray implementation of IMAGE converts arguments to full

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


    varargin = cellfun(@full,varargin,'UniformOutput',false);
    segwarn(SegArray,'converting SegArray to full for IMAGE function');
    try
        [varargout{1:nargout}] = image(varargin{:});
    catch me
        throwAsCaller(me);
    end
end
