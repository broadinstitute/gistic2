function p=add_patches(v,ctable,pw,pl,ps,lw)



% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

ax=axis;
axis([ax(1:2) ax(3)-ps-pl ax(4)]);

if ~exist('lw','var')
  lw=0.1;
end

for i=1:length(v)
  if ~isnan(v(i))
    if v(i)==0
      col=[ 1 1 1];
      alp=0;
    else
      ci=clip_to_range(round(v(i)),[1 size(ctable,1)]);
      col=ctable(ci,:);
      alp=1;
    end
    if ischar(col)
        p(i)=text(i,ax(3)-ps-pl/2,deblank(col),'HorizontalAlignment','center',...
            'VerticalAlignment','baseline','FontSize',lw); 
    else
        p(i)=patch([i-pw/2 i+pw/2 i+pw/2 i-pw/2 i-pw/2],[ax(3)-ps ax(3)-ps ax(3)-ps-pl ax(3)-ps-pl ax(3)-ps], ...
                col);
        set(p(i),'FaceAlpha',alp);
        set(p(i),'LineStyle','-');
        if lw==0
            set(p(i),'LineStyle','none','LineWidth',eps);
        else
            set(p(i),'LineWidth',lw);
        end
    end
  else
    p(i)=NaN;
  end
end

