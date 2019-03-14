function [arm_medians median_vector centromere_snps] = find_med_chrarm(CL,cyto)
% FIND_MED_CHRARMS - find median copy levels of chromosome arms

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

  armnames = unique(cellfun(@char,regexp({cyto.name},'^[0-9]+[p-q]+','match'),'UniformOutput',false));
  armnames = armnames(2:end);
  
  % make chromosome arms structure
  chrarms = struct('name',{},'start',{},'end',{},'chrn',{},'length',{});
  for i=1:length(armnames)
    idx = strmatch(armnames(i),{cyto.name});
    band.name = char(armnames(i));
    band.start = cyto(idx(1)).start+1;
    band.end = cyto(idx(end)).end;
    band.chrn=cyto(idx(end)).chrn;
    band.length = band.end-band.start+1;
    chrarms(i) = band;
  end
  % sort the arms by chromosome
  [~,order]=sort([chrarms.chrn]);
  chrarms = chrarms(order);
  
  % find the centromeres
  Ncyto_chrn = ceil(length(chrarms)/2);
  centromeres = zeros(Ncyto_chrn);
  for i=1:Ncyto_chrn
    centromeres(i) = chrarms(2*i-1).end;
  end

  Ndata_chrn = max(CL.chrn);
  arm_medians = zeros(2*Ndata_chrn,size(CL.dat,2));
  centromere_snps = NaN(1,Ndata_chrn);
  disp('Finding median of arms on chromosome:')
  for ch=1:Ndata_chrn
    disp(ch)

    ch_snps = SegArray(CL.chrn == ch);
    above_cent = SegArray(CL.pos>centromeres(ch));
    p_snps = ch_snps & ~above_cent;
    q_snps = ch_snps & above_cent;
    if any(p_snps)
        arm_medians(2*ch-1,:) = median(CL.dat(p_snps,:),1);
        centromere_snps(ch) = find(p_snps,1,'last');
    else
        arm_medians(2*ch-1,:) = NaN;
    end
    if any(q_snps)
        arm_medians(2*ch,:) = median(CL.dat(q_snps,:),1);
    else
        arm_medians(2*ch,:) = NaN;
    end
    
  end
  median_vector = arm_medians';
  median_vector = median_vector(:)';
  
    
