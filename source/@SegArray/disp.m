function disp(S)
% overloaded SegArray DISP method

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

if ~isempty(S.vals)
    gapsize = 4;
    termsize = get(0,'CommandWindowSize');
    termwide = termsize(1);

    % map row position to text position
    rowmap = cumsum(any(S.bpts,2));
    % find brkrows and gaps
    [I,J] = find(S.bpts);
    brkrows = unique([I(:);size(S.bpts,1)]);
    NR = length(brkrows);
    gaps = diff(brkrows) > 1;
    % gap-marked space
    gaptext = repmat(' ',NR,4);
    gaptext(gaps,:) = '_';
    % initialize output array as row index
    rowindex = [num2str(brkrows),repmat(':',NR,1)];
    charmat = rowindex;
    % initial output array = row index
    startcol = 1;
    % loop over columns
    for c = 1:size(S.bpts,2)
        % create a text column of data values
        valtext = num2str(S.vals(J==c));
        coltext = repmat([repmat(' ',1,size(valtext,2)-1),'"'],NR,1);
        coltext(rowmap(I(J==c)),:) = valtext;
        coltext = [gaptext coltext];
        % if no more room in window, emit columns
        if size(charmat,2) + size(coltext,2) >= termwide
            % output as many columns as will fit
            fprintf('\nColumns %d through %d\n',startcol,c-1);
            disp(charmat);
            charmat = rowindex;
            startcol = c;
        end
        charmat = [charmat coltext]; 
    end
    % output remaining columns
    if startcol > 1
        fprintf('\n Columns %d through %d\n',startcol,c);
    end
    disp(charmat);
end


