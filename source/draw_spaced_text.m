function draw_spaced_text(t,a,ymin,ymax,x1,x2,x3,f,extend_dense,lineparams,varargin)
% t - text
% a - actual position
% ymin, ymax - range of all figure
% x1, x2, x3 - posotion of arm and elbow and text (0.1, 0.3, 0.32)

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


axis([0 1 ymin ymax]);
set(gca,'YDir','reverse');
n=length(t);
a=as_column(a);

if exist('varargin','var')
    th=text(x1*ones(n,1),a,t,varargin{:});
else
    th=text(x1*ones(n,1),a,t);
end

h=get(th,'Extent');
if(iscell(h))
  h=cat(1,h{:});
else
  h=cat(1,h);
end
w=max(h(:,4));
s=w*f;

if exist('extend_dense','var') && extend_dense
    d=dist(a,[]);
    d=d+diag(nan(n,1));
    idx1=find(min(d,[],1)<w);
end

y=find_text_pos(a,w,s,ymin,ymax);
idx=find(abs(y-a)<s*0.1);
if ~isempty(y)
    y(idx)=a(idx);
end

for i=1:length(th)
    p=get(th(i),'Position');
    p(2)=y(i);
    set(th(i),'Position',p);
end

if exist('extend_dense','var') && extend_dense
    idx=union(find(a~=y),idx1);
else
    idx=1:length(y);
end
   
for j=1:length(idx)
    i=idx(j);
    p=get(th(i),'Position');
    p(1)=x3;
    set(th(i),'Position',p);    
    if exist('lineparams','var') && ~isempty(lineparams)
        lh(i)=line([0 x1 x2],[a(i) a(i) y(i)],'Color','k','LineWidth',0.5,lineparams{:});
    else
        lh(i)=line([0 x1 x2],[a(i) a(i) y(i)],'Color','k','LineWidth',0.5);
    end
end

return
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% draw_spaced_text({'aaaa','bb','CCCC','DDD'},[10 50 53 55],1,100,0.2,0.5,0.55,1,1);
% draw_spaced_text({'aaaa','bb','CCCC','DDD'},[10 50 53 55],1,100,0.1,0.3,0.32,0.5,0,[],'FontSize',9);
% Flipped view
% clf; draw_spaced_text({'aaaa','bb','CCCC','DDD'},[10 50 53 55],1,100,0.1,0.3,0.32,0.5,0,[],'FontSize',9,'HorizontalAlignment','Right'); set(gca,'XDir','reverse')

