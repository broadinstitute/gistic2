function out = vertcatfill(varargin)
%out = horzcatfill(varargin)

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


sizes = cellfun(@size,varargin,repmat({2},1,length(varargin)));
endsizes = max(sizes).*ones(1,length(sizes));

if ~isequal(sizes,endsizes)

    varargin = cellfun(@emptyfill,varargin,mat2cell(endsizes,1,ones(1,length(endsizes))),'UniformOutput',0);

end

out = horzcat(varargin);
    

function fill = emptyfill(x,l)

if isnumeric(x(1)) || ischar(x(1))
    fill = horzcat(x,cast(zeros(1,l-size(x,2)),class(x(1))));
elseif iscell(x(1)) 
    fill = horzcat(x,repmat({[]},1,l-size(x,2)));
elseif isstruct(x(1))
    emptystruct = x(1);
    for fl = fieldnames(emptystruct)'
        emptystruct.(char(fl)) = '';
    end
    fill = horzcat(x,repmat(emptystruct,1,l-size(x,2)));
end

    