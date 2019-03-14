function draw_colorbar(horiz,posvec,textvec,fs)


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

if strcmp(horiz(1:4),'hori')
  if (0)
    imagesc((1:64));  
    set(gca,'XTick',posvec,'XTickLabel',textvec,'FontSize',fs,'YTick',[]);
    set(gca,'TickLength',[0 0]);
  else
    pos=get(gca,'Position');
    pos2=[pos(1) pos(2)+pos(4)/2 pos(3) pos(4)/2];
    set(gca,'Position',pos2);
    imagesc((1:64)); 
    for i=1:length(posvec)
        text(posvec(i),2,textvec{i},'FontSize',fs,'HorizontalAlignment','center');
    end 
    set(gca,'XTick',posvec,'XTickLabel',[]);
    set(gca,'YTick',[],'XTickLabel',[]);    
    set(gca,'TickLength',[0 0]);
%    set(gca,'TickLength',[0.005 0.0125]);
  end
else
    imagesc((1:64)');
    set(gca,'YDir','normal');
%    axis([0 2 0.5 64.5]);
%    for i=1:length(posvec)
%        text(posvec(i),0.5,textvec{i},'FontSize',fs);
%    end 
    set(gca,'YTick',posvec,'YTickLabel',textvec,'FontSize',fs,'XTick',[]);
%    set(gca,'TickLength',[0.005 0.0125]);
    set(gca,'TickLength',[0 0]);
end
