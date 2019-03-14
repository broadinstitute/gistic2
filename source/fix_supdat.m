function D=fix_supdat(D)


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

if iscell(D.supacc)
  is_cell=1;
else
  D.supacc=cellstr(D.supacc);
  D.supdesc=cellstr(D.supdesc);
  is_cell=0;
end

for i=1:size(D.supdat,1)
  if isempty(find(D.supacc{i}==':'))
    u=unique(D.supdat(i,:));
    u=sort(u(~isnan(u)));
    if ~isempty(setdiff(u,[0 1]))
      supacc=[ D.supacc{i} ': ' ];
      supdesc=[ D.supdesc{i} ': ' ];
      v=D.supdat(i,:);
      for j=1:length(u)
        if j<length(u)
          st='/';
        else
          st='';
        end
        v(D.supdat(i,:)==u(j))=j;
        supacc=[supacc num2str(j) '-[' num2str(u(j)) ']' st];
        supdesc=[supdesc num2str(j) '-[' num2str(u(j)) ']' st];
      end
      D.supacc{i}=supacc;
      D.supdesc{i}=supdesc;
      D.supdat(i,:)=v;
    end
  end
end

if ~is_cell
  D.supacc=strvcat(D.supacc);
  D.supdesc=strvcat(D.supdesc);
end
