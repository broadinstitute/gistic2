function [fract chrarms]= normalize_by_arm_length(D,Q,cyto,norm_type,ref_length,chrarms,skip1)
  
%% Function takes in a segmented data structure and its associated
%ziggurat matrix Q and two arrays, one containing
%the armlengths of each chr_arm (sorted by chr arm) and the start
%positions of each chr_arm.  It returns an array of length size(Q,1)
%containing the fractional length of each segment normalized by the
%length of the chromosome arm on which it sits.  If a segment crosses
%the centromere, then the function divides by the length of the arm on
%which the majority of the segment sits.
%  
%norm_type: 1 = by number of snps (default); 2 = by length in bp
%ref_length: 1 = max % of either chr arm (i.e. max value is 1); 2 = sum
%of % of chromosome arms (i.e. max value is 2)
%skip1: if chr arms not interrogated, max length for segments on chr = 1
%when skip1 = 0 (default); segments of length 1 are made = NaN when skip1 = 1
%---
% $Id$
% $Date: 2014-01-31 15:30:18 -0500 (Fri, 31 Jan 2014) $
% $LastChangedBy: schum $
% $Rev$

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


  if ~exist('norm_type','var')
    norm_type = 1;
  end
 
  if ~exist('ref_length','var')
    ref_length = 1;
  end
  
  if ~exist('skip1','var')
    skip1=0;
  end
  
  if ~exist('chrarms','var') || isempty(chrarms)   
    armnames = unique(cellfun(@char,regexp({cyto.name},'^[0-9 X Y]+[p-q]+','match'),'UniformOutput',false)); 
%    armnames = unique(cellfun(@char,regexp({cyto.name},'^[0-9]+[p-q]+','match'),'UniformOutput',false)); 
%    armnames = armnames(2:end);
    chrarms = {};
    skip1chr = [];
    count=0;
    for i=1:length(armnames)
      idx = strmatch(armnames(i),{cyto.name});
      band.name = char(armnames(i));
      band.start = cyto(idx(1)).start+1;
      band.end = cyto(idx(end)).end;
      band.chrn=cyto(idx(end)).chrn;
      band.length = band.end-band.start+1;
      bnd_snps = find_snps(D,band.chrn,band.start,band.end,0);
      band.snp_start = min(bnd_snps);
      band.snp_end = max(bnd_snps);
      band.snp_length = length(bnd_snps);
      if band.snp_length==0;
        band.snp_start = -1;
        band.snp_end = -1;
        skip1chr=unique([skip1chr band.chrn]);
      end
      if(length(find(D.chrn==band.chrn))>0)
        count=count+1;
        chrarms{1}(count) = band;
      end
    end
    clear count
    
    chr=cat(1,chrarms{1}.chrn);  
    [spos,sposi]=sort(chr);
    chrarms{1} = chrarms{1}(sposi);
    
  end
  
  armlengths = [chrarms{1}.length];
  armstart = [chrarms{1}.start];
  armlengths_by_snp = [chrarms{1}.snp_length];
  armstart_by_snp = [chrarms{1}.snp_start];
  
  isq = find(D.pos(Q(:,2))>=armstart(2*Q(:,1))');
  spans_cent = intersect(find(D.pos(Q(:,2))<armstart(2*Q(:,1))'), ...
                         find(D.pos(Q(:,3))>=armstart(2*Q(:,1))'));
  isp = setdiff(1:size(Q,1),union(isq,spans_cent));  
  
  QQ = Q;
  for i=1:length(spans_cent)
    if mod(i,1000) == 0
      disp(i)
    end
    switch ref_length
      case {1}
       if (armstart_by_snp(2*Q(spans_cent(i),1)'-Q(spans_cent(i),2))/armlengths_by_snp(2*Q(spans_cent(i),1)'-1) > ...
             (Q(spans_cent(i),3)-armstart_by_snp(2*Q(spans_cent(i),1)'))/armlengths_by_snp(2*Q(spans_cent(i),1)'))
         isp(length(isp)+1) = spans_cent(i);
         Q(spans_cent(i),3) = armstart_by_snp(2*Q(spans_cent(i),1)');
       else
         isq(length(isq)+1) = spans_cent(i);
         Q(spans_cent(i),2) = armstart_by_snp(2*Q(spans_cent(i),1)');
       end
     case {2}
        isp(length(isp)+1) = spans_cent(i);
        Q(spans_cent(i),3) = armstart_by_snp(2*Q(spans_cent(i),1)')-1;
        QQ(spans_cent(i),2) = armstart_by_snp(2*Q(spans_cent(i),1)');
    end
  end
  
  norms = zeros(1,size(Q,1));
  
  switch norm_type
   case {1}
    norms(isq) = armlengths_by_snp(2*Q(isq,1)');
    norms(isp) = armlengths_by_snp((2*Q(isp,1)-1)');
    lengths = Q(:,3) - Q(:,2)+1;
    fract = lengths./norms';
    if ref_length > 1
      if ~isempty(spans_cent)
        norms(spans_cent) = armlengths_by_snp(2*QQ(spans_cent,1)');
        lengths(spans_cent) = QQ(spans_cent,3)-QQ(spans_cent,2)+1;
        fract(spans_cent) = fract(spans_cent)+lengths(spans_cent)./norms(spans_cent)';
      end 
    end
   case {2}
    norms(isq) = armlengths(2*Q(isq,1)');
    norms(isp) = armlengths((2*Q(isp,1)-1)');
    lengths = (D.pos(Q(:,3))-D.pos(Q(:,2)))+1;
    fract = lengths./norms';
    if ref_length > 1
      if ~isempty(spans_cent)
        norms(spans_cent) = armlengths(2*QQ(spans_cent,1)');
        lengths(spans_cent) = QQ(spans_cent,3)-QQ(spans_cent,2)+1;
        fract(spans_cent) = fract(spans_cent)+lengths(spans_cent)./norms(spans_cent)';
      end 
    end
  end

  if(skip1==1)
    for i = 1:size(Q,1)
      if(ismember(Q(i,1),skip1chr) & fract(i)==1)
        fract(i) = NaN;
      end
    end
  end
  
