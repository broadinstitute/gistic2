function [termA, termAstr] = get_labels(stuff)
%
% ---
% $Id$
% $Date: 2014-01-31 15:29:30 -0500 (Fri, 31 Jan 2014) $
% $LastChangedBy: schum $
% $Rev$

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


if isempty(stuff)
    stuff = '0.25';
end

counter = 1;
for k = 1:size(stuff,2)
    if ~isspace(stuff(k))
        stuff2(counter) = stuff(k);
        counter = counter +1;
    end
end

if stuff2(end) == ','
    stuff2(end) = '';
end

found = find(stuff2==',');

if isempty(found)
    termA(1) = str2num(stuff2(1:size(stuff2,2)));
	termAstr{1} = stuff2(1:size(stuff,2));
else
    for k = 1:size(found,2)
        if k == 1
            termA(k) = str2num(stuff2(1:found(k)-1));
            termAstr{k} = stuff2(1:found(k)-1);
        elseif k > 1 && k <= size(found,2)
            termA(k) = str2num(stuff2(found(k-1)+1:found(k)-1));
            termAstr{k} = stuff2(found(k-1)+1:found(k)-1);
        end
    end
    termA(size(found,2)+1) = str2num(stuff2(found(k)+1:end));
    termAstr{size(found,2)+1} = stuff2(found(k)+1:end);
end

for k = 1:size(termAstr,2)
    found2 = find(termAstr{k} == '^');
    if ~isempty(found2)
        termAstr{k} = [termAstr{k}(1:found2) '{' termAstr{k}(found2+1:end) '}'];
    end
end
