function th=draw_names_box(nms,horiz,fs,vec,stagger,baseline,textparams)


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

if isempty(nms)
  th=NaN;
  return
end

if iscell(nms)
  nms=strvcat(nms);
end

if (nargin<4) || isempty(vec)
  vec=ones(size(nms,1),1);
end

if nargin<5
  stagger=0;
end
cs=0.5*(cumsum([0; vec(1:(end-1))])+cumsum(vec))+0.5;

if exist('textparams','var')
  parpos=find(strcmp(textparams,'vec_is_pos'));
  if ~isempty(parpos)
    cs=vec;
    textparams(parpos)=[];
  end
end

cs2=nan(size(cs));
if exist('textparams','var')
  parpos=find(strcmp(textparams,'vec2'));
  if ~isempty(parpos)
    cs2=textparams(parpos+1);
    textparams(parpos:(parpos+1))=[];
  end
elseif stagger
end

set(gca,'YDir','reverse');
if strcmp(lower(horiz(1:3)),'hor')
  axis([0.5 sum(vec)+0.5 0 1]);
  if ~exist('baseline','var') || isempty(baseline)
    baseline=0.25;
  end
  for i=1:size(nms,1)
    if exist('textparams','var')
      th(i)=text(cs(i),baseline+stagger*mod(i,2),deblank(nms(i,:)),...
                 'Rotation',90,'VerticalAlignment','middle','HorizontalAlignment',...
                 'center','FontSize',fs,'Interpreter','none',textparams{:});
    else
      th(i)=text(cs(i),baseline+stagger*mod(i,2),deblank(nms(i,:)),...
                 'Rotation',90,'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',fs,'Interpreter','none');
    end
  end
else %vertical
  axis([0 1 0.5 sum(vec)+0.5]);
  if ~exist('baseline','var') || isempty(baseline)
    baseline=1;
  end
  for i=1:size(nms,1)
    if exist('textparams','var')
      th(i)=text(baseline+stagger*mod(i,2),cs(i),deblank(nms(i,:)),'HorizontalAlignment','Right',...
                 'VerticalAlignment','Middle','FontSize',fs, ...
                 'Interpreter','none',textparams{:}); 
    else
      th(i)=text(baseline+stagger*mod(i,2),cs(i),deblank(nms(i,:)),'HorizontalAlignment','Right',...
                 'VerticalAlignment','Middle','FontSize',fs, ...
                 'Interpreter','none'); 
    end
%    e=get(th(i),'Extent');
%    lh(i)=line([e(1)+e(3) 0.03],[ cs(i) cs(i)],'Color',[0 0 0]);
%    rh(i)=rectangle('position',[0.03 cs(i)-vec(i)*0.5 0.01 vec(i)],'curvature',[0.3 0.3]);
  end
end
