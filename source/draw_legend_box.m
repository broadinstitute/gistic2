function draw_legend_box(D,ncols,fst,fsw,title_indent,box_width,varargin)


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

if ~exist('title_indent','var')
  title_indent=0;
end

if ~exist('box_width','var')
  box_width=0.5;
end

nl=0;
tw=1.5;
for i=1:size(D.supdat,1)
  [st,sn]=break_sup_names(D.supacc(i,:));
  if isempty(sn)
      continue;
  end   
  nl=nl+tw; % title
  seen=unique(D.supdat(i,:));
  seen(isnan(seen))=[];
  nl=nl+ceil(length(seen)/ncols);
end   

axis([0 ncols 0 nl+0.5]);
set(gca,'YDir','reverse');

nl=0;
for i=1:size(D.supdat,1)
  [st,sn]=break_sup_names(D.supacc(i,:));
  if isempty(sn)
      continue;
  end   
  nl=nl+tw; % title
  text(title_indent,nl-tw/3,st,'FontSize',fst);
  seen=unique(D.supdat(i,:));
  seen(isnan(seen))=[];
  if ~isempty(seen) && seen(1)==0
      add_seen=1;
  else
      add_seen=0;
  end   
  for j=1:length(seen)
      c=mod(j-1,ncols);
      nl=nl+(c==0);
      if seen(j)==0
        patch([c+0.05 c+box_width c+box_width c+0.05 c+0.05],...
              [nl-0.1 nl-0.1 nl-0.9 nl-0.9 nl-0.1],[ 1 1 1],'FaceAlpha',0);
      else
        if isempty(D.supmark(i).colormap)
          hold on;
          ph=plot(c+0.25,nl-0.5,'x');
          %'MarkerSize',D.supmark(i).marker_size(seen(j)),...
          set(ph,'Marker',D.supmark(i).marker{seen(j)},...
                 'MarkerSize',fsw*2.5, ...
                 'LineWidth',fsw*0.1,...
                 'MarkerFaceColor',[ 0 0 0],'MarkerEdgeColor',[ 0 0 ...
                              0]);
          axis off;
        else
          patch([c+0.05 c+box_width c+box_width c+0.05 c+0.05],...
                [nl-0.1 nl-0.1 nl-0.9 nl-0.9 nl-0.1],...
                D.supmark(i).colormap(seen(j),:));
        end 
      end
      if (seen(j)+add_seen)>length(sn)
        txt='EMPTY';
      else
        txt=sn{seen(j)+add_seen};
      end
      th=text(c+box_width+0.05,nl-0.5,txt,'FontSize',fsw,'HorizontalAlignment','left','Interpreter','none');
%      if ~isempty(varargin)
%        set(th,varargin{:});
%      end
  end   
end   
