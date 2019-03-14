function [res,gr]=display_D(Dord,gdend,sdend,disp_p,in_pos)


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

if ~exist('gdend','var')
  if isfield(Dord,'gdend')
    gdend=Dord.gdend;
  else
    gdend=[];
  end
end

if ~exist('sdend','var')
  if isfield(Dord,'sdend')
    sdend=Dord.sdend;
  else
    sdend=[];
  end
end

if exist('disp_p','var')
  if ischar(disp_p)
    disp_p=get_disp_defaults(disp_p);
  elseif iscell(disp_p)
    if ischar(disp_p{1})
      disp_p1=get_disp_defaults(disp_p{1});
    else
      disp_p1=disp_p{1};
    end
    disp_p1=add_struct(disp_p1,disp_p{2});
    disp_p=disp_p1;
  end
else
  disp_p=get_disp_defaults;
end

if exist('in_pos','var') && ~isempty(in_pos)
  gr=make_subplotgrid(disp_p.x.sizes,disp_p.y.sizes,disp_p.x.gaps, ...
                      disp_p.y.gaps,disp_p.x.border,disp_p.y.border,in_pos);
else
  gr=make_subplotgrid(disp_p.x.sizes,disp_p.y.sizes,disp_p.x.gaps, ...
                      disp_p.y.gaps,disp_p.x.border,disp_p.y.border);  
end

for i=1:length(disp_p.items)
  [h,gr]=subplotgrid(gr,disp_p.items{i}{1},disp_p.items{i}{2});
  if length(disp_p.items{i})>3
    res{i}=display_d_elem(gr,disp_p.items{i}{3},Dord,gdend,sdend, ...
                   disp_p.items{i}{4:end});
  else
    res{i}=display_d_elem(gr,disp_p.items{i}{3},Dord,gdend,sdend);
  end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
%       SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is a "subfunctioned" and reduced version of the the display_D_elem
% function. It was created to decouple GISTIC source code dependence from
% preprocess_D and about 40 called functions.
%
function res=display_d_elem(gr,elem_type,Dord,gdend,sdend,varargin)
% display element in current axes

res=[];
switch elem_type
 case 'dataorig' %! keep
  % sample the copy number vertically
  maxypix = 10000; % maximum number of pixels
  samplint = ceil(size(Dord.dat,1)/maxypix);
  % draw heatmap
  imagesc(Dord.dat(1:samplint:end,:));
  if ~isempty(varargin)
    caxis(varargin{1});
  end
  % no axis ticks
  set(gca,'XTick',[]);
  set(gca,'YTick',[]);
  axis on;
  box on;
  bluepink;
  if ~isempty(varargin)
    if length(varargin)>1
      colormap(varargin{2});
    end
  end
  res=caxis;
 
 case 'ssupdat' %! keep
  if isfield(Dord,'supdat')
    if ~isfield(Dord,'supmark')
      c=read_colorscheme();
      Dord=add_supmark(Dord,c);
    end
    if nnz(Dord.supdat<0)
      disp('Fixing supdat');
      Dord=fix_supdat(Dord);
    end
    axis([0.5 size(Dord.supdat,2)+0.5 0 0.000001]);
    for i=1:size(Dord.supdat,1)
      if ~isfield(Dord.supmark(i),'patchwidth')
        pw=0.8;
      else
        pw=Dord.supmark(i).patchwidth;
      end
      if isfield(Dord.supmark(i),'linewidth')
        add_patches(Dord.supdat(i,:),Dord.supmark(i).colormap,pw, ...
                    Dord.supmark(i).height,0.5*(i>1),Dord.supmark(i).linewidth);
      else
        add_patches(Dord.supdat(i,:),Dord.supmark(i).colormap,pw, ...
                    Dord.supmark(i).height,0.5*(i>1));
      end
    end
    noticks;
  end
  
 case 'ssupacc' %! keep
  if isfield(Dord,'supdat')
    if ~isfield(Dord,'supmark')
      c=read_colorscheme();
      Dord=add_supmark(Dord,c);
    end
    hvec=cat(1,Dord.supmark(:).height)+0.5;
    hvec(1)=hvec(1)-0.25;
    hvec(end)=hvec(end)-0.25;
    [suptitle,supnames]=break_sup_names(Dord.supacc);
    if length(varargin)>2
      res=draw_names_box(strvcat(suptitle),varargin{1},varargin{2},hvec,varargin{3:end});
    else
      res=draw_names_box(strvcat(suptitle),varargin{1},varargin{2},hvec);
    end
  end
  noticks;
  
 case 'ssuplegend' %! keep
  if isfield(Dord,'supdat')
    if ~isfield(Dord,'supmark')
      c=read_colorscheme();
      Dord=add_supmark(Dord,c);
    end    
    draw_legend_box(Dord,varargin{:});
  end


 case 'colorbar' %! keep
  switch length(varargin)
   case 0
    draw_colorbar('hori');
   case 4
    draw_colorbar(varargin{1},varargin{2},varargin{3},varargin{4})
   case 6 % was 5 
    if strcmp(varargin{2},'for')
      r=get_subplotgrid_data_range(gr,varargin{3},varargin{4});
      draw_colorbar(varargin{1},[1+(63/varargin{5})*(0:varargin{5})],cellstr(num2str((r(1):(r(2)-r(1))/varargin{5}:r(2))'))',varargin{6})      
    else
      error('unknown colorbar param');
    end
   otherwise
      error('wrong number of parameters for colorbar');    
  end   
%  noticks;
  axis on
  box on 

  case 'chrn' %! keep
   if isfield(Dord,'chrn')
      [u,ui,uj]=unique_keepord(Dord.chrn);
      rl=runlength(uj);
      sz=rl(:,2)-rl(:,1)+1;
      res=draw_names_box(num2chromosome(unique_keepord(Dord.chrn)),varargin{1},varargin{2},sz,varargin{3:end});
   end 
  case 'chrcyto' %! keep
   if exist('varargin','var') && length(varargin)>0
     if isfield(Dord,'cyto_stain')
       % take care of NaNs
       tmp=Dord.cyto_stain;
       tmp(isnan(tmp))=0;
       tmp=tmp/100;
       image(repmat(1-tmp,[1 1 3]));
       noticks;
       axis on
       box on
     end
   elseif isfield(Dord,'chrn')
        image(repmat(mod(Dord.chrn,2)*0.9,[1 1 3]));
        if isfield(Dord,'armn')
          ax=axis;
          armpos=find(diff(Dord.armn)>0);
          for i=1:length(armpos)
            res(i)=line(ax(1:2),[armpos(i) armpos(i)]+0.5,'LineStyle','-','Color',[0.6 0.6 0.6],'LineWidth',1);
          end
        end
        noticks;            
   end
  
end
