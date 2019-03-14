function display(S)
% overloaded SegArray DISPLAY method

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

    
    try
        % test internal representation before displaying
        validate(S,inputname(1));
    
        if ~isempty(S.vals)
            % snippet below taken from 'help display'
            if isequal(get(0,'FormatSpacing'),'compact')
                disp([inputname(1) ' =']);
                disp(S);
            else
                disp(' ');
                disp([inputname(1) ' =']);
                disp(' ');
                disp(S);
            end
        else
            fprintf(1,'Empty SegArray: %d-by-%d\n\n',size(S,1),size(S,2));
        end
    catch me
        % exception thrown by SegArray validation
        disp(['Invalid SegArray: ' me.message]);
    end
