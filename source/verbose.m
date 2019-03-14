function verbose(str,level,varargin)
% VERBOSE display a message according to the global verbose level
%    VERBOSE(STR,LEVEL,VARARGIN) displays STR if the LEVEL is greater or equal
%    to the current global VERBOSE_LEVEL.  VARARGIN can be used to pass
%    arguments to formatting commands in STR.  
% 
%    USe these standard levels:
%           10 for a general messege 
%           20 for more details
%           30 for most extensive information
%
%   Examples:
%
%       verbose('General message',10) displays 'General message' 
%       only if the VERBOSE_LEVEL is 10 or above
%
%       verbose('Today''s date is %s',10,datestr(now,'yymmdd')) displays
%       today's date iff the VERBOSE_LEVEL is 10 or above.
%
%
%
%    See also SET_VERBOSE_LEVEL
%
%    History:
%
%        26 Sept 07:  Changed disp to fprintf to allow for escape
%        formatting.  (Jen Dobson, jdobson@broad.mit.edu)
%
%        5 Nov 07:  Added regexprep for % characters so that remainder of
%        line would not be commented out.
%
%        27 Nov 07:  Added option to pass formating strings.
%
%---
% $Id$
% $Date$
% $LastChangedBy$
% $Rev$

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.



global VERBOSE_LEVEL
global VERBOSE_FILE

if nargin==1
    level=1;
end


if isempty(varargin)

str = char(regexprep(cellstr(str),'%','%%'));
    if ~isempty(VERBOSE_LEVEL) & (level<=VERBOSE_LEVEL)
        fprintf(1,[str repmat('\n',size(str,1),1)]');  %escape the % to prevent line from commenting
        if ~isempty(VERBOSE_FILE)
            if ~exist(VERBOSE_FILE,'file')
                fid = fopen(VERBOSE_FILE,'w');
            else
            fid = fopen(VERBOSE_FILE,'a');
            end
            
            fprintf(fid,str,varargin{:});  %escape the % to prevent line from commenting
            fprintf(fid,'\n');
            fclose(fid);
        end
    end

else

    if ~isempty(VERBOSE_LEVEL) & (level<=VERBOSE_LEVEL)
        fprintf(1,str,varargin{:});  %escape the % to prevent line from commenting
        fprintf(1,'\n')
        if ~isempty(VERBOSE_FILE)
            fid = fopen(VERBOSE_FILE,'a');
            fprintf(fid,str,varargin{:});  %escape the % to prevent line from commenting
            fprintf(fid,'\n');
            fclose(fid);
        end
    end
end
