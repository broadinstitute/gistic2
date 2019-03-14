function gr=make_subplotgrid(vx,vy,spx,spy,bx,by,in_pos)


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

set(gcf,'renderer','zbuffer');
% set(gcf,'renderer','none');

if ~exist('in_pos','var')
  in_pos=[0 0 1 1]; % lrx,lry,w,h
  clf;
end

asp=(sum(vx)/in_pos(3))/(sum(vy)/in_pos(4));

%if asp>(8/10.5) % 0.25 from each side
if asp>1
  set(gcf,'Position',[50 50 1000 1000/asp]);
%  set(gcf,'PaperSize',[11 8.5]);
  set(gcf,'PaperPosition',[0.25 1.5 8 8/asp]);
%  set(gcf,'PaperSize',[22 17]);
%  set(gcf,'PaperPosition',[0.25 1.5 16.5 16.5/asp]); 
else
  set(gcf,'Position',[50 50 800*asp 800]);
%  set(gcf,'PaperSize',[11 8.5]);
  set(gcf,'PaperPosition',[0.25 1.5 8*asp 8]); 
%  set(gcf,'PaperSize',[22 17]);  
%  set(gcf,'PaperPosition',[0.25 1.5 20*asp 20]);
end

vx=vx./sum(vx)*(1-bx)*in_pos(3);
vy=vy./sum(vy)*(1-by)*in_pos(4);
spx=spx./sum(spx)*bx*in_pos(3);
spy=spy./sum(spy)*by*in_pos(4);

nx=length(vx);
ny=length(vy);

gr=cell(ny,nx);
sy=in_pos(2)+in_pos(4); % was 1
for i=1:ny
  sy=sy-spy(i);
  sx=in_pos(1); % was 0
  for j=1:nx
    sx=sx+spx(j);
    gr{i,j}.position=[ sx sy-vy(i) vx(j) vy(i) ]; 
    gr{i,j}.handle=NaN;
%    gr{i,j}.handle=subplot('position',gr{i,j}.position);
%    axis off
%    box off
    sx=sx+vx(j);
  end
  sy=sy-vy(i);
end

