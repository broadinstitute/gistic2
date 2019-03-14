function rl=replace_rl(rl,rli,val)


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

if size(rl,2)<4
  rl=[rl zeros(size(rl,1),1)];
end

if iscell(rl)
  for i=1:length(rl)
    if ~isempty(rli{i})
      rl{i}=replace_rl(rl{i},rli{i},val{i});
    end
  end
else
  rl(rli,3)=val;
  
  if (0)
    while ~isempty(rli)
      % connect with prev
      if rli(1)>1 && rl(rli(1)-1,3)==val && (rl(rli(1)-1,4)==rl(rli(1),4))
        en=rl(rli(1),2);
        rl=rl(setdiff(1:size(rl,1),rli(1)),:);
        rli=rli-1;
        rl(rli(1),2)=en;
      end
      
      % connect with next
      if rli(1)<size(rl,1) && rl(rli(1)+1,3)==val && (rl(rli(1),4)==rl(rli(1)+1,4))
        en=rl(rli(1)+1,2);
        rl=rl(setdiff(1:size(rl,1),rli(1)+1),:);
        rl(rli(1),2)=en;
        rli(2:end)=rli(2:end)-1;
      end
      rli(1)=[];
    end
  end
  
  rl_of_rl=runlength(rl(:,3),rl(:,4));
  rl_of_rl(:,4)=rl(rl_of_rl(:,1),4);
  new_rl=zeros(size(rl_of_rl,1),4);
  for j=1:size(new_rl,1)
    new_rl(j,1)=rl(rl_of_rl(j,1),1);
    new_rl(j,2)=rl(rl_of_rl(j,2),2);
    new_rl(j,3:4)=rl_of_rl(j,3:4);
  end
  if any(new_rl(:,2)-new_rl(:,1)<0) || any(new_rl(:,1)==0) || any(new_rl(:,2)==0)
    error('cannot be');
  end
  rl=new_rl;
end
