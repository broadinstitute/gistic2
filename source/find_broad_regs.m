function regs=find_broad_regs(D,cyto,regs,broad_type)


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

if ~isfield(D,'chrn')
  D=add_chrn(D);
end

if ischar(cyto)
  cyto=read_cytoband_data(cyto);
end

if ~isfield(D,'cyton')
  D=add_cyto(D,cyto);
end

if ~exist('broad_type','var')
  broad_type=struct('method','region_size','p_arm',0.5);
elseif isnumeric(broad_type)
  broad_type=struct('method','region_size','p_arm',broad_type);
end

switch broad_type.method
 case 'region_size'
  prev_st=0;
  prev_en=0;
  % identify broad regions
  for k=1:length(regs)
    for i=1:length(regs{k})
      if (regs{k}(i).st==prev_st) && (regs{k}(i).en==prev_en)
        regs{k}(i).broad=0;
        regs{k}(i).broad_ratios=-1;
        continue;
      end
      c=D.chrn(regs{k}(i).peak_st);
      c_pos=find(D.chrn==c);
      arm_length=histc(D.armn(c_pos),1:2);
      reg_arms=histc(D.armn(regs{k}(i).st:regs{k}(i).en),1:2);
      rat=as_column(reg_arms)./(arm_length+eps);
      if any(rat>=broad_type.p_arm)
        regs{k}(i).broad=1;
      else
        regs{k}(i).broad=0;
      end
      regs{k}(i).broad_ratios=rat;
      prev_st=regs{k}(i).st;
      prev_en=regs{k}(i).en;
    end
  end
 case 'peeloff'
  prev_st=0;
  prev_en=0;
  % identify broad regions
  for k=1:length(regs)
    for i=1:length(regs{k})
      if (regs{k}(i).st==prev_st) && (regs{k}(i).en==prev_en)
        regs{k}(i).broad=0;
        regs{k}(i).broad_ratios=-1;
        continue;
      end
      c=D.chrn(regs{k}(i).peak_st);
      c_pos=find(D.chrn==c);
      arm_length=histc(D.armn(c_pos),1:2);
      for j=1:length(regs{k}(i).segments)
        if ~isempty(regs{k}(i).segments{j})
          reg_arms=histc(D.armn(regs{k}(i).segments{j}.st:regs{k}(i).segments{j}.en),1:2);
          rat=as_column(reg_arms)./(arm_length+eps);
          if any(rat>=broad_type.p_arm)
            regs{k}(i).segments{j}.broad=1;
          else
            regs{k}(i).segments{j}.broad=0;
          end
          
        end
      end
      segs=cat(1,regs{k}(i).segments{:});
      b=cat(1,segs.broad);
      regs{k}(i).broad_ratios=sum(b)/length(b);
      if regs{k}(i).broad_ratios>broad_type.freq
        regs{k}(i).broad=1;
      else
        regs{k}(i).broad=0;
      end       
      prev_st=regs{k}(i).st;
      prev_en=regs{k}(i).en;
    end
  end      
  
 case 'scorecutoff'
  % identify broad regions
  ads=broad_type.ads;
  for k=1:length(regs)
    for i=1:length(regs{k})
      %      if (regs{k}(i).st==prev_st) && (regs{k}(i).en==prev_en)
      %        regs{k}(i).broad=0;
      %        regs{k}(i).broad_ratios=-1;
      %        continue;
      %      end
      c=D.chrn(regs{k}(i).peak_st);
      c_pos=find(D.chrn==c);
      c_st=min(c_pos);
      arm_length=histc(D.armn(c_pos),1:2);
	  % handle peaks that violate platau assumption (thak maximal value as
	  % score)
      score_at_peak=max(ads{k}(regs{k}(i).peak_st:regs{k}(i).peak_en),[],1);  
      
      if broad_type.score_thresh(k)<broad_type.score_thresh_focal(k)
        error('cannot be');
      end
      if score_at_peak<broad_type.score_thresh(k)
        error('cannot be such a low peak');
      end
      % check for broad
      if score_at_peak<=(broad_type.score_thresh(k)+broad_type.score_thresh_focal(k))
        sc_cutoff=score_at_peak*(broad_type.score_thresh_focal(k))/(broad_type.score_thresh_focal(k)+ ...
                                                          broad_type.score_thresh(k));
        if isfield(broad_type,'min_cutoff')
          sc_cutoff=max(sc_cutoff,broad_type.min_cutoff(k));
        end
        x=ads{k}(c_pos)>=sc_cutoff;
        
        xrl=runlength(x);
        rli=find(xrl(:,1)<=regs{k}(i).peak-c_st+1 & xrl(:,2)>=regs{k}(i).peak-c_st+1);
        reg_arms=histc(D.armn(c_pos(xrl(rli,1):xrl(rli,2))),1:2);
        rat=as_column(reg_arms)./(arm_length+eps);
        if any(rat>=broad_type.p_arm)
          regs{k}(i).broad=1;
          regs{k}(i).broad_ratios=rat;
          regs{k}(i).broad_st=xrl(rli,1)+c_st-1;
          regs{k}(i).broad_en=xrl(rli,2)+c_st-1;
          regs{k}(i).focal=0;
          regs{k}(i).focal_ratios=rat;          
        elseif all(rat<broad_type.p_arm)
          regs{k}(i).focal=1;
          regs{k}(i).focal_ratios=rat;          
          regs{k}(i).broad=0;
          regs{k}(i).broad_ratios=rat;
        end
      else   
        reg_arms=histc(D.armn(regs{k}(i).st:regs{k}(i).en),1:2);
        rat=as_column(reg_arms)./(arm_length+eps);
        if any(rat>=broad_type.p_arm)
          regs{k}(i).broad=1;
          regs{k}(i).broad_ratios=rat;
          regs{k}(i).broad_st=regs{k}(i).st;
          regs{k}(i).broad_en=regs{k}(i).en;
        else
          regs{k}(i).broad=0;
          regs{k}(i).broad_ratios=rat; 
        end
        
        % check for focal
        x=ads{k}(c_pos)>=score_at_peak-broad_type.score_thresh_focal(k);
        xrl=runlength(x);
        rli=find(xrl(:,1)<=regs{k}(i).peak-c_st+1 & xrl(:,2)>=regs{k}(i).peak-c_st+1);
        reg_arms=histc(D.armn(c_pos(xrl(rli,1):xrl(rli,2))),1:2);
        rat=as_column(reg_arms)./(arm_length+eps);
        if all(rat<broad_type.p_arm)
          regs{k}(i).focal=1;
          regs{k}(i).focal_ratios=rat;          
        else
          regs{k}(i).focal=0;
          regs{k}(i).focal_ratios=rat;
        end
      end
    end
  end      
  
% case 
end

