function sp = plot_snp_score(base_dir,ext,D,q,ads,q_thresh,use_loglog,write_pdf,~,...
    cyto,qv_scale,qv_ticks,score_ticks,xsz,flip_score_pv,no_stripes,regs,add_broad_bars,...
    add_peak_names,genepattern,make_png,flip_dels,gene_name_side,fname)
%PLOT_SNP_SCORE Plot q values for amplifications and deletions
%
%   SP = PLOT_SNP_SCORE(BASE_DIR, FNAME_SUFFIX, D, Q, ADS, Q_THRESH, ...
%                       USE_LOGLOG, WRITE_PDF, unused, CYTO, QV_SCALE,...
%                       QV_TICKS, SCORE_TICKS, XSZ, FLIP_SCORE_PV,...
%                       NO_STRIPES, REGS, ADD_BROAD_BARS, ADD_PEAK_NAMES,...
%                       GENEPATTERN, MAKE_PNG, FLIP_DELS, GENE_NAME_SIDE)
%
% required arguments
%   BASE_DIR - directory for image file outputs
%   FNAME_SUFFIX - suffix to "amplification" or "deletion"
%   D - structure containing copy number data
%   Q - a two-element amp/del cell array of q-value vectors for
% amplifications and deletions respectively
%   ADS - a two-element amp/del cell array of amplification and deletion
% G-scores
%   Q_THRESH - the significance threshold for Q values
%   USE_LOGLOG controls how the qvalues/scores are scaled:
%           0 = plot log10(qvalue)
%           1 = plot log10(-log10(qvalue)+1)
%           2 = plot log10(score) %!!! broken
%           3 = plot score
%
% optional parameters
%   WRITE_PDF - boolean, set to write plot figures to the output directory
%   UNUSED - place-holder for a formerly important parameter
%   CYTO - cytoband information
%   QV_SCALE - two-element amp/del vector specifying the maximum value of 
% the q-value scale (USE_LOGLOG = 1 or 2)
%   QV_TICKS - two-element amp/del array of structs with fields 'val' for
% q-values where labels should be places and 'txt' for the label text 
%   SCORE_TICKS - two-element amp/del array of structs with fields 'val' for
% G-scores where labels should be places and 'txt' for the label text
%   XSZ - width of plot? (default = 1)
%   FLIP_SCORE_PV - boolean, set to 1 to put q-values along top(default=0)
%   NO_STRIPES - - boolean, set to 1 to supress chromosome bar image
%   REGS - two-element amp/del cell array of struct arrays representing peak 
%  regions
%   ADD_BROAD_BARS - boolean set to mark broad regions on the plot
% (currently disabled).
%   ADD_PEAK_NAMES - TODOC remaining pareameters %!!!
%   GENEPATTERN -
%   MAKE_PNG -
%   FLIP_DELS -
%   GENE_NAME_SIDE -
%

% GISTIC software version 2.0
% Copyright (c) 2011, 2016 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


%   1 July 08: Added 'close all' at end of each loop on length(q) --
%   jdobson@broad.mit.edu
%
% $Id$
% $Date$
% $LastChangedBy$
% $Rev$

% add_wide_peaks,
%{
    for i=1:max(C.chrn)
        chrnpos(i)=round(mean(find(C.chrn==i)));
    end
%}

if ~exist('flip_dels','var') || isempty(flip_dels)
  flip_dels=0;
end

if ~exist('gene_name_side','var') || isempty(gene_name_side)
  gene_name_side=[1 1-2*(flip_dels>0)]; % 1 is right, -1 is left
end

if ~exist('xsz','var') || isempty(xsz)
    xsz=[1 1];
end

if ~exist('write_pdf','var') || isempty(write_pdf)
    write_pdf = true;
end

if ~exist('make_png','var') || isempty(make_png)
    make_png = false;
end

if ~exist('flip_score_pv','var') || isempty(flip_score_pv)
    flip_score_pv=[0 0];
end
if length(flip_score_pv)==1
    flip_score_pv=[ flip_score_pv flip_score_pv];
end

if ~exist('no_stripes','var') || isempty(no_stripes)
    no_stripes=0;
end

if ~exist('add_broad_bars','var') || isempty(add_broad_bars)
    add_broad_bars=0;
end

if ~exist('fname','var')
    fname = '';
end

% add arm information if we have it
if exist('cyto','var')
    c=combine_cyto(cyto,'arm');
    D=rmfield_if_exists(D,{'cyton','armn','cyto_stain'});
    D=add_cyto(D,c);
end

%% loop over amplifications(1), then deletions (2)

% name (extension) for image files
names={'amp_qplot','del_qplot'};
% line color for q/score
line_cols='rb';
% interpolation data
sort_q = cell(1,2);
ord_q = cell(1,2);
sort_ads = cell(1,2);
ord_ads = cell(1,2);

for k=1:2
  % initialize interpolation data
  [sort_q{k},order] = sort(q{k});
  ord_ads{k} = ads{k}(order);
  
  % add extrapolated point for end of scale
  if use_loglog==1
      diffpt = find(sort_q{k}~=sort_q{k}(1),1,'first');
      if ~isempty(diffpt)
          lq1 = log10(sort_q{k}(1));
          lq2 = log10(sort_q{k}(diffpt));
          lqx = log10(qv_ticks(k).vals(end));
          adx = ord_ads{k}(1)+(ord_ads{k}(diffpt)-ord_ads{k}(1)).*(lqx-lq1)/(lq2-lq1);
          ord_ads{k} = [adx;ord_ads{k}];
          sort_q{k} = [10^lqx;sort_q{k}];
      end
  end
  %!!! TODO find out why ads(q) is not monotonically nonincreasing
  [sort_ads{k},order] = sort(ord_ads{k});
  ord_q{k} = sort_q{k}(order);
  
  name = [base_dir fname names{k} ext];
  
  % start figure
  figure; clf;
  
  if gene_name_side(k)==1
    gr=make_subplotgrid([0.2 0.06 xsz(k) 0.4 ],[0.05 2 0.05],[1 1 0.3 0 1],[1 0 0 1],0.1,0.1);
  else
    gr=make_subplotgrid([0.2 0.06 0.4 xsz(k) ],[0.05 2 0.05],[1 1 0.3 0 1],[1 0 0 1],0.1,0.1);
  end
  
  % draw q-values
  if gene_name_side(k)==1
    sp{k}=subplotgrid(gr,2,3);
  else
    sp{k}=subplotgrid(gr,2,4);
  end    
  %  sp{k}=subplot(10,1,1:9);
  if ~exist('qv_scale','var') || isempty(qv_scale)
    for kk=1:2
      qv_scale(kk)=0.9*min(q{kk});
    end
  end
  if ~exist('qv_ticks','var')
    qv_ticks=[];
  end
  if ~exist('score_ticks','var')
    score_ticks=[];
  end
  
  %% establish scaling for q-value/score
  ymin=0;
  minq = 1e-320; % minimum loggable q-value in matlab's floating point representation
  if use_loglog==3 % plot score and not q-value
    y = ads{k};
    ymax = max(max(y),max([score_ticks(k).vals]));
    ymax = ymax*1.05;
  elseif use_loglog==2 % log of score
                       % find minimal score
    ymin=log10(max(ads{k}(q{k}>0.99)));
    y=log10(ads{k});
    % if use_loglog==2 then qv scale is actually score scale
    ymax=log10(qv_scale(k));
    if ~isempty(qv_ticks)
      qv_ticks(k).vals=log10(qv_to_score(qv_ticks(k).vals,sort_q,ord_ads,k));
    end
    % CONTINUE HERE 
  elseif use_loglog==1
    % log (negative) log q scale
    qscale = @(q) log10(-log10(max(q,minq))+1);
    y = qscale(q{k});
    ymax = qscale(qv_ticks(k).vals(end));
%!    ymax = qscale(qv_scale(k));
  else %! use_loglog==0
    % use negative log q scale
    qscale = @(q) -log10(max(q,minq));
    y = qscale(q{k});
    ymax = qscale(qv_scale(k));
  end
  %-> y = scaled q-value or score for plotting
  %-> ymax = upper range for score
  %-> qscale is scale function for q-value
  
  xcoor=1:length(q{k});
  ycoor=y;
  for ci=1:max(D.chrn)
    in_chr=find(D.chrn==ci);
    xcoor=[xcoor min(in_chr)-0.3 max(in_chr)+0.3];
    ycoor=[ycoor; ymin-1; ymin-1];
  end
  [sx,sxi]=sort(xcoor);
  xcoor=sx;
  ycoor=ycoor(sxi);
  plot(xcoor,ycoor,'Color',line_cols(k)); hold on
  y_gt_ymax=find(ycoor>ymax);
  % draw lines at edge where score/q is off scale
  if ~isempty(y_gt_ymax)>0
    plot(xcoor(y_gt_ymax),ymax*ones(length(y_gt_ymax),1),'m.'); % [line_cols(k) '.']); 
  end
  ax=axis;
  axis([0.5 length(q{k})+0.5 ymin ymax]); %min(ax(4),30)]);
  ax=axis;
  set(gca,'XTick',[]);
  set(gca,'CameraUpVector',[-1 0 0]);
  set(gca,'YTick',[]);

  if flip_dels && k==2
    camera_position=get(gca,'CameraPosition');
    camera_position(end)=-camera_position(end);
    set(gca,'CameraPosition',camera_position);
  end
  
  %% create shaded chromosome blocks
  for ci=2:2:max(D.chrn)
    inchr=find(D.chrn==ci);
    if ~isempty(inchr)
      mnc=min(inchr);
      mxc=max(inchr);
      ph=patch([ mnc-0.5 mxc+0.5 mxc+0.5 mnc-0.5 mnc-0.5],...
               [ ax(3) ax(3) ax(4) ax(4) ax(3)],[ 0.9 0.9 0.9]);
      set(ph,'FaceAlpha',0.9,'EdgeColor','none');
    end
  end
  box on
  %% indicate centromere locations
  if isfield(D,'cyton')
    for j=grep('q',{c.name},1)
      b=find(D.cyton==j,1,'first');
      if ~isempty(b)
        %           disp(c(j).name);
        %           disp(b);
        lh=line([b-0.5 b-0.5],[ax(3) ax(4)],'LineStyle',':');
        set(lh,'Color',[0.3 0.3 0.3]);
      end
    end
  end
  % redraw plot (why?)
  plot(xcoor,ycoor,'Color',line_cols(k));
  set(gca,'CameraUpVector',[-1 0 0]);
  
  ax=axis;
  axis([0.5 length(q{k})+0.5 ax(3:4)]);
  ax=axis;
  line([ ax(1) ax(2) ax(2) ax(1) ax(1)],[ax(3) ax(3) ax(4) ax(4) ax(3)],'Color',[0 0 0]);
  
  %% draw score and q-value scales
  
  geneL = length(q{k});     % genomic length (in markers)
  txt_opts = {'HorizontalAlignment','center','Rotation',0};
  darkgreen=[91 186 71]/255; % color for significance threshold

  x1a = [0.5 geneL*0.005];  % upper scale tick extents
  x1b = 0.5-geneL*0.001;    % upper scale text location
  x1b_align = 'bottom';     % upper scale alignment
   
  x2a = [0.5+geneL geneL*(1-0.005)]; % lower scale text extents
  x2b = 0.5+geneL*1.005;     % lower scale text location
  x2b_align = 'top';         % lower scale alignment

  % if exchanging scales, exchange alignments
  if flip_score_pv(k)
    [x1a,x2a]=exchange_vars(x1a,x2a);
    [x1b,x2b]=exchange_vars(x1b,x2b);
    [x1b_align,x2b_align]=exchange_vars(x1b_align,x2b_align);
  end
  
  if use_loglog==3 % plot score, not q-value
      % render the independant score scale
      if ~isempty(ads) && ~isempty(score_ticks)
        for i=1:length(score_ticks(k).vals)
          val = score_ticks(k).vals(i);
          if val <= ymax;
            line(x1a,[val val],'Color','k');
            text(x1b,val,score_ticks(k).txt{i},txt_opts{:},'VerticalAlignment',x1b_align);
          end
        end
      end
    % the q-value scale depends on the score scale
    if ~isempty(qv_ticks)
      for i=1:length(qv_ticks(k).vals)
        tmp=find(q{k}<=qv_ticks(k).vals(i));
        if ~isempty(tmp)
          val=min(y(tmp));
          line(x2a,[val val],'Color','k');
          text(x2b,val,qv_ticks(k).txt{i},txt_opts{:},'VerticalAlignment',x2b_align);
        end
      end
    end
  else % use_loglog ~= 3
    % render the independant q-value scale
    if ~isempty(qv_ticks)
      for i=1:length(qv_ticks(k).vals)
        val = qscale(qv_ticks(k).vals(i));
        line(x1a,[val val],'Color','k');
        text(x1b,val,qv_ticks(k).txt{i},txt_opts{:},'VerticalAlignment',x1b_align);
      end
    end
    % render the score scale (depends on the q-value score)
    if ~isempty(ads) && ~isempty(score_ticks)
        for i=1:length(score_ticks(k).vals)
          if use_loglog==2
            line(x2a,log10(score_ticks(k).vals(i))*ones(1,2),'Color','k');
            text(x2b,log10(score_ticks(k).vals(i)),score_ticks(k).txt{i},txt_opts{:},'VerticalAlignment',x2b_align);
          else % use_loglog is 0 or 1
            % calculate the position by converting ads to q
            % and then scaling to the plot 
            if score_ticks(k).vals(i) < sort_ads{k}(end)
              val = qscale(score_to_qv(score_ticks(k).vals(i),sort_ads,ord_q,k));
              line(x2a,[val val],'Color','k');
              text(x2b,val,score_ticks(k).txt{i},txt_opts{:},'VerticalAlignment',x2b_align);
            end
          end
        end
    end
  end
  
  %% significance threshold
  axs(1)=gca;
  if q_thresh > 0
    if use_loglog==3
      tmp=find(q{k}<=q_thresh);
      if ~isempty(tmp)
        val_green=min(y(tmp));
        line([1 geneL],[val_green val_green],'Color',darkgreen);
      end
    elseif use_loglog==2
      line([1 geneL],[1 1]*log10(qv_to_score(q_thresh,sort_q,ord_ads,k)),'Color',darkgreen);      
    else % use_loglog is 0 or 1
      line([1 geneL],[1 1]*qscale(q_thresh),'Color',darkgreen);
    end
  end
    
  %% chromosome numbers
  subplotgrid(gr,2,1);
  [u,ui,uj] = lunique(D.chrn);
  
  %Write Chromosome numbers on plot
  draw_names_box(num2chromosome(u),'ver',10,histc(uj,1:length(ui)),-0.5,1);
  box off
  axis off
  
  %% draw stripes
  if ~no_stripes
    subplotgrid(gr,2,2);  
    % subplot(10,1,10);
    colormap gray
    imagesc(mod(D.chrn,2));
    axs(2)=gca;
    set(gca,'XTick',[],'YTick',[]);
  end

  %% add broad and focal bars
  if exist('regs','var') && ~isempty(regs) && add_broad_bars
    subplotgrid(gr,2,2);
    %    axis([0 1 0.5 length(q{k})+0.5]);
    axis([0 2 0.5 length(q{k})+0.5]);
    set(gca,'YDir','reverse');
    pcols=[ 1 0 0; 0 0 1];
    
    for i=1:length(regs{k})
      if regs{k}(i).broad
        if isfield(regs{k}(i),'broad_st')
          bph=patch([ 0.1 0.9 0.9 0.1 0.1],[regs{k}(i).broad_st regs{k}(i).broad_st regs{k}(i).broad_en ...
                              regs{k}(i).broad_en ...
                              regs{k}(i).broad_st],pcols(k,:));
        else
          bph=patch([ 0.1 0.9 0.9 0.1 0.1],[regs{k}(i).st regs{k}(i).st regs{k}(i).en regs{k}(i).en ...
                              regs{k}(i).st],pcols(k,:));
        end          
        set(bph,'EdgeColor',[1 1 1],'LineWidth',0.25);
      end
      if regs{k}(i).focal 
        %        bph=patch([ 1.2 1.9 1.9 1.2 1.2],[regs{k}(i).peak_wide_st regs{k}(i).peak_wide_st ...
        %                            regs{k}(i).peak_wide_en regs{k}(i).peak_wide_en ...
        %                            regs{k}(i).peak_wide_st],pcols(k,:));
        %        set(bph,'EdgeColor',pcols(k,:),'LineWidth',0.25);
        line([1.2 1.9],regs{k}(i).peak*ones(1,2),'LineWidth',0.5,'Color',pcols(k,:));
      end      
    end
    text(0.5,length(q{k})*1.005,'Broad','Rotation',90,'FontSize',8,'HorizontalAlignment','Right','VerticalAlignment','middle');
    text(1.75,length(q{k})*1.005,'Focal','Rotation',90,'FontSize',8,'HorizontalAlignment','Right','VerticalAlignment','middle');
    axs(2)=gca;
    set(gca,'XTick',[],'YTick',[]);
  end
  
  %% add some kind of wide peak markings (forced off)
  if 0 && exist('regs','var') && ~isempty(regs) && add_wide_peaks
    subplotgrid(gr,2,4);
    axis([0 1 0.5 length(q{k})+0.5]);
    set(gca,'YDir','reverse');
    pcols=[ 1 0 0; 0 0 1];
    n=length(q{k});
    for i=1:length(regs{k})
      if add_wide_peaks==1 %wide_peak
        y11=regs{k}(i).peak_wide_st-0.5;
        y21=regs{k}(i).peak_wide_en+0.5;
        y12=0.5*(y11+y21)-n*0.001;
        y22=0.5*(y11+y21)+n*0.001;
        y13=y11-n*0.005;
        y23=y21+n*0.005;
        bph=patch([ 0.1 0.5 0.5 0.9 0.9 0.5 0.5 0.1],[y11 y13 y12 y12 y22 y22 y23 y21],pcols(k,:));
      elseif add_wide_peaks==2 % peak
        y11=regs{k}(i).peak-0.5;
        y21=regs{k}(i).peak+0.5;
        y12=0.5*(y11+y21)-n*0.001;
        y22=0.5*(y11+y21)+n*0.001;
        y13=y11-n*0.005;
        y23=y21+n*0.005;
        %        disp([0.5*(y12+y22)-0.5*(y13+y23)])
        bph=patch([ 0.1 0.5 0.5 0.9 0.9 0.5 0.5 0.1],[y11 y13 y12 y12 y22 y22 y23 y21],pcols(k,:));
      end
      
      %       set(bph,'EdgeColor',[1 1 1],'LineWidth',0.001); % 'LineStyle','none',
      set(bph,'LineStyle','none');
      set(bph,'Clipping','off');
    end
    axs(2)=gca;
    set(gca,'XTick',[],'YTick',[]);
  end
  
  %% label the peaks
  if exist('regs','var') && ~isempty(regs) && exist('add_peak_names','var') && ~isempty(add_peak_names)
    if gene_name_side(k)==1
      subplotgrid(gr,2,4);
    else
      subplotgrid(gr,2,3);
    end      
    axis([0 1 0.5 length(q{k})+0.5]);
    set(gca,'YDir','reverse');
    if genepattern
      warning('off','optim:quadprog:SwitchToMedScale')
    end
    
    if ~isempty(add_peak_names{k})
      tmp=cat(1,regs{k}.peak);
      if isfield(regs{k},'added_manually')
        tmp(cat(1,regs{k}.added_manually)==1)=[]; % support adding peaks at the end that are not in the
                                                        % all_lesions_file
      end
      if gene_name_side(k)==1
        draw_spaced_text(add_peak_names{k},tmp,1,length(q{k}),0.1,0.3,0.32,0.5,0,[],...
                         'FontSize',9);
      else
        draw_spaced_text(add_peak_names{k},tmp,1,length(q{k}),0.1,0.3,0.32,0.5,0,[],...
                         'FontSize',9,'HorizontalAlignment','Right'); 
        set(gca,'XDir','reverse');
      end        
    end
    warning('on','optim:quadprog:SwitchToMedScale');
  end
  
%!    set(gca,'XTick',chrnpos,'XTickLabel',C.chr(chrnpos));
    
  %% write image file(s)

  file_formats = [repmat({'pdf'},write_pdf),...
                  repmat({'png'},write_pdf & make_png),...
                  repmat({'fig'},~genepattern)];
  save_current_figure(name,file_formats);
%{
  if exist('write_pdf','var') && write_pdf
    %  print_D(name,{{'fig'}});
    %  print_pdf(name);
    
    % use vector graphics drivers to get MS quality output
    rend = get(gcf,'renderer'); % save renderer for later restore
    if str2double(regexprep(version,'\.[0-9]+\.[0-9]+ .+$','')) >= 8.4
        % new (matlab R2014b on) vector graphics output method
        set(gcf,'renderer','painters');
        cf = get(gcf,'Number');
    else
        % old (pre-matlab R2014b) vector graphics output method
        set(gcf,'renderer','none');
        cf = num2str(gcf);
    end
    %!print(['-f' cf],'-depsc',[name,'.eps']);
    print(['-f' cf],'-dpdf',[name,'.pdf']);
    set(gcf,'renderer',rend); % restore renderer
    if make_png
      saveas(gcf,[name '.png'],'png');
    end
    %   print('-depsc',[name '.eps']); 
  end
  if ~genepattern
    saveas(gcf,[name '.fig'],'fig');
  end
%}
  
  close all

end

%% subfunction - map q-values to score
function y=qv_to_score(x,q,ads,k)
  y = zeros(1,length(x));
  for i=1:length(x)
    y(i) = interp_pwl(q{k},ads{k},x(i));
  end
%% subfunction - map score to q-value  
function x = score_to_qv(y,ads,q,k)
  x = zeros(1,length(y));
  for i=1:length(x)
      x(i) = 10.^interp_pwl(ads{k},log10(q{k}),y(1));
  end
