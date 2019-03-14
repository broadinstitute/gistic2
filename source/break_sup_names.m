function [suptitle,supnames]=break_sup_names(nms)


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

if ischar(nms)
  nms=cellstr(nms);
end

for i=1:size(nms,1)
    posacc=find(nms{i}==':'); 
    if isempty(posacc)
        suptitle{i}=deblank(nms{i});
        supnames{i}={};
        continue;
    end
    suptitle{i}=nms{i}(1:(posacc(1)-1));
    
    sacc=posacc(1)+2;
    pacc=find([ nms{i} '/']=='/');
    assert(~isempty(pacc));
    
    for j=1:length(pacc)
        cur=nms{i}(sacc:(pacc(j)-1));
        dashpos=find(cur=='-');
        if ~isempty(dashpos)
            cur=cur((dashpos(1)+1):end);
        end
        curnames{j}=cur;     
        sacc=pacc(j)+1;
    end 
    supnames{i}=curnames;
end 

if size(nms,1)==1
    suptitle=suptitle{1};
    supnames=supnames{1};
end 

