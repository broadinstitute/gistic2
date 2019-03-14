function save_current_figure(fname,exts)
%Save matlab figures with various formats
%
% save_figure(FNAME,EXTS)
% FNAME is the required base name for the figure (including path)
% EXTS is one of 'fig', 'png', 'pdf' or 'eps' (default 'fig') for
%   extension/file type

% GISTIC software version 2.0
% Copyright (c) 2011, 2016 Gad Getz, Rameen Beroukhim, Craig Mermel, Jen Dobson,
% Steve Schumacher, Nico Stransky, Mike Lawrence, Gordon Saksena
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

if ~exist('exts','var') || isempty(exts)
    exts = {'fig'};
end

% change only one argument into a list
if ~iscellstr(exts)
    exts = {exts};
end

% fig = matlab format
if ismember('fig',exts)
    saveas(gcf,[fname,'.fig'],'fig');
end

% png = WWW standard raster format
if ismember('png',exts)
    saveas(gcf,[fname,'.png'],'png');
end

% pdf = portable document format (vector) 
if ismember('pdf',exts)
    [cf,oldrend] = vector_render();
    print(['-f' cf],'-dpdf',[fname,'.pdf']);
    set(gcf,'renderer',oldrend); % restore renderer
end

% eps = encapsulated post script vector format
if ismember('eps',exts)
    [cf,oldrend] = vector_render();
    print(['-f' cf],'-depsc',[fname,'.eps']);
    set(gcf,'renderer',oldrend); % restore renderer
end

function [cf,oldrend] = vector_render()
    % test version for new matlab graphics
    oldrend = get(gcf,'renderer'); % save renderer for later restore
    new_graphix = str2double(regexprep(version,'\.[0-9]+\.[0-9]+ .+$','')) >= 8.4;
    if new_graphix
        set(gcf,'renderer','painters');
        cf = num2str(get(gcf,'Number'));
    else
        set(gcf,'renderer','none');
        cf = num2str(gcf);
    end
