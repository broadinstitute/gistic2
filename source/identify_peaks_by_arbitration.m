function [regs seg_assignment] = identify_peaks_by_arbitration(Z,Qs,ads,q,...
    score_thresh,gene_gistic_rg,gene_gistic_ds,gene_gistic_res,do_arbitration,arm_peeloff)
%IDENTIFY_PEAKS_BY_ARBITRATION - find significant regions on the genome
%
%   [REGS SEG_ASSIGNMENT] = IDENTIFY_PEAKS_BY_ARBITRATION(Z,QS,ADS,Q,...
%         SCORE_THRESH,GENE_GISTIC_RG,GENE_GISTIC_DS,GENE_GISTIC_RES,...
%         DO_ARBITRATION)
%
% input:
%   Z  - a "D" structure with chrn, pos, and dat fields. Z.dat is a two
%        element cell array: Z.dat{1} are amplification scores and Z.dat{2}
%        are deletion scores.
%   QS - focal ziggurat segments for
%   ADS - two-element cell array of amp/del scores for each marker
%   Q - amp/del cell of vectors of marker-based q-values 
%   SCORE_THRESH - score at which (marker-based) q-value becomes
%       significant, 1x2 vector for [amps,dels]
%   GENE_GISTIC_RG - reference genome with fields added by gene gistic:
%       'snps' are the marker indices for snps representing the gene, 
%       'sc' is the gene score (maximum of scores of markers in the gene), 
%       'p' and 'q' are significance statistics of the gene score.
%   GENE_GISTIC_DS - a cell array of empirical score distributions for each
%       gene size (in snp units). Each distribution is a normalized vector
%       of probabilities for binned scores. The width of each score bin is
%       set by the GENE_GISTIC_RES input.
%   GENE_GISTIC_RES - resolution of the score bins used in the 
%       distributions of gene_gistic_ds
%   DO_ARBITRATION - optional flag set (default) to use arbitrated peeloff
%
% output:
%   REGS - two-element amp/del cell array of struct arrays representing peak 
%       regions that are the result of this stage of the analysis. These
%       contains all the information about that peak except for the wide
%       peak regions which are added by calculate_wide_peaks_for_regs.
%       Each struct in the region arrays contains the following fields
%         'chrn' chromosome number
%         'samples'  1xN vector containing 1 if that sample has an alteration
%             covering the peak, 0 if not 
%         'peak_st' & 'peak_en' the the peak limits (the narrow peak or mcr),
%             in marker space 
%         'peak'  probe closest to the middle of the peak limits
%         'qv' original q-value assigned to the peak marker by GISTIC
%         'st' & 'en' limits of the region that exceeds qv_thresh around
%             the peak, in marker space 
%         'orig_score'  original score assigned to the peak marker by GISTIC
%         'iter_score'  score assigned to the peak marker at the
%             beginning of arbitrated peel-off for that peak 
%         'score' the score assigned to the peak marker at the end of
%             arbitrated peel-off 
%         'ads' the score for that chromosome at the end of arbitrated
%             peel-off (used by RegBounder to calculate wide peaks) 
%         'resid_qv'  the q-value associated with score
%         'robust_reg_st' & 'robust_reg_en'  the starting and ending point
%             of the boundaries set by score-G_t, where G_t is the score
%             associated with a q-value of qv_thresh.  This is the region
%             within which Regbounder will look for peak limits.
%
%   SEG_ASSIGNMENT - two element amp/del cell array of assignment matrices.
%       Each assignment matrix shows how the score from each event (row) is
%       distributed across peak regions (columns). The row-wise sum of an
%       assignment matrix would equal the event score input; the
%       column-wise sum would equal the peak score output.
%

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


%! TODO: has do_arbitration = false been tested?

  if ~exist('do_arbitration','var') || isempty(do_arbitration)
    do_arbitration = 1;
  end
  
  if ~exist('arm_peeloff','var') || isempty(arm_peeloff)
    arm_peeloff = false;
  end
  
  if ~exist('gene_gistic_rg','var') || isempty(gene_gistic_rg)
    do_gene_gistic = 0;
  else
    do_gene_gistic = 1;
    if ~exist('gene_gistic_ds','var') || isempty(gene_gistic_ds)
      error(['Must supply gene gistic background distributions to use gene ' ...
             'gistic']);
    end
  end
  
  use_segarray = isa(Z.dat{1},'SegArray');
  
  regs = cell(1,2);
  seg_assignment = cell(1,2);
  
  %% loop across SCNA type: amps, then dels
  for k=1:2
    if k==1
      Q = Qs.amp;
    else
      Q = Qs.del;
    end
    
    % identify event arms for arm peeloff
    if arm_peeloff
        Q(:,10) = NaN;
        % find all unique arms represented in the markers
        arm_rl = runlength(Z.armn,Z.chrn);
        arm_rl(:,3) = 2 * Z.chrn(arm_rl(:,1)) + arm_rl(:,3) - 2;
        % identify the arm to which each event belongs in column 10
        for j=1:size(arm_rl,1)
            Q(Q(:,2)>=arm_rl(j,1) & Q(:,3)<=arm_rl(j,2),10) = arm_rl(j,3);
        end
    end
    
    orig_Q = Q;
    orig_ads = ads{k};
    nsegs = size(Q,1);
    Q(:,11) = (1:nsegs)';
    
    n=length(regs{k});
    
    if k==2 && do_gene_gistic
      gene_gistic = 1;
      gene_scores = [gene_gistic_rg.gene_scores];
      gene_q = [gene_gistic_rg.q];
      [sorted_ads,sadsi] = sort(gene_scores);
      sorted_q = log(gene_q(sadsi));
    else
      gene_gistic = 0;
      [sorted_ads,sadsi] = sort(ads{k});
      sorted_q = log(q{k}(sadsi));
    end
    
    %% loop over chromosomes
    for ch=1:max(Z.chrn)
      loop = 0;
      disp(['Finding peaks on chromosome ' num2str(ch)])
      % Find snps within that chromosome
      if use_segarray
        in_chr = SegArray(Z.chrn)==ch; 
        chr_zero = find(in_chr,1,'first')-1; % snp offset of current chromosome
        chr_max = find(in_chr,1,'last');
       else
        in_chr = find(Z.chrn==ch);
        chr_zero = min(in_chr)-1; % Marks first snp not on current chromosome
        chr_max = max(in_chr); % Marks last snp on current chromosome
      end
      regs_on_chr = [];
      x=Z.dat{k}(in_chr,:); % x contains only data from present
                            % chromosome
      y=(x>0);
      ns=size(x,2); % ns = # samples
      
      if gene_gistic
        orig_ads = nanmean(x,2);
        genes_in_chr = find([gene_gistic_rg.chrn] == ch);
        % gene gistic score already stored in reference genome
        sc = [gene_gistic_rg(genes_in_chr).gene_scores];
        orig_sc = sc;
        [mx gene_mi] = max(sc);
        top_genes = find(sc == mx);
        % rg = snps covering top_gene; mid snp of top gene
        mi = gene_gistic_rg(genes_in_chr(gene_mi)).midsnp-chr_zero;
        % get marker offsets within chromosome for top genes
        rg = [];
        for i=1:length(top_genes)
          rg = [rg; gene_gistic_rg(genes_in_chr(top_genes(i))).snps-chr_zero];
        end
        rg = unique(rg);
      else  
        % snp score is mean across samples
        sc = nanmean(x,2); % take row-wise mean of each data
        orig_sc = sc;
        [mx mi] = max(sc); % find maximum score and max snp
        rg=find(sc==mx); % find entire peak region where sc = mx

        % pick larger of any noncontiguous regions
        [rg,mi] = fix_rg(rg,mi);
      end

      if mx >= score_thresh(k)
        loop = 1;
      end
      
      %% main peak identification loop
      while loop
        n=n+1; % n = current reg number
                
        % Find all segs covering this peak
        all_segs = find_snp_segments(orig_Q,rg+chr_zero,arm_peeloff);
        %disjoint_segs = find_snp_segments(Q,rg+chr_zero);
        
        % remove the peak from the data
        [Znew Qnew Qrm] = ziggurat_peeloff(Z.dat{k},Q,rg+chr_zero,[],arm_peeloff);
        disjoint_segs = Qrm(:,11); % indices of removed events
        
        % Create seg_assignment column for new reg
        seg_assignment{k}(:,n) = zeros(1,nsegs);
        % For each disjoint reg, assign all of height to this reg
        seg_assignment{k}(disjoint_segs,n) = Qrm(:,4);
        
        % Identify conflicting seg (over peak but already peeled off)
        conflict_segs = setdiff(all_segs,disjoint_segs);
                
        %% Resolve conflicts
        if ~isempty(conflict_segs) && do_arbitration
          disp('Resolving conflicts...')
          % calculate total score of unconflicting events
          no_conflict_segs = setdiff(1:nsegs,conflict_segs);
          disjoint_weight = sum(seg_assignment{k}(no_conflict_segs,:),1)/size(Z.dat{k},2);
          
          dist_weights = zeros(length(conflict_segs),n);
          new_sc = sc;
          % for each conflicting event...
          for i=1:length(conflict_segs)
            % identify peaks that the conflicting event belongs to
            conflicting_regs = find(seg_assignment{k}(conflict_segs(i),:)>0);
            conflicting_regs = union(conflicting_regs,n);
            % weight these peaks by their sizes
            % (measured with unconflicting events) 
            reg_wgts = disjoint_weight(conflicting_regs);
            reg_wgts = reg_wgts / sum(reg_wgts);
            % distribute each event's score to peaks according to weight
            dist_weights(i,conflicting_regs) = reg_wgts * orig_Q(conflict_segs(i),4);
            cur_amt_added = dist_weights(i,n);
            % adjust scores
            cur_seg = orig_Q(conflict_segs(i),:);
            if gene_gistic
              if ~use_segarray
                Z.dat{k}(cur_seg(2):cur_seg(3),cur_seg(5)) = ...
                        Z.dat{k}(cur_seg(2):cur_seg(3),cur_seg(5)) + cur_amt_added;
              end
            else
              new_sc(cur_seg(2)-chr_zero:cur_seg(3)-chr_zero) = ...
                  new_sc(cur_seg(2)-chr_zero:cur_seg(3)-chr_zero) + cur_amt_added/size(Z.dat{k},2);
            end
          end
          
          % if gene gistic, recalculate gene scores
          if gene_gistic
            if use_segarray
              % adjust scores
              Z.dat{k} = addSegments(Z.dat{k},orig_Q(conflict_segs,2),...
                            orig_Q(conflict_segs,3),orig_Q(conflict_segs,5),...
                            dist_weights(:,n));
            end
            %
            conflict_snps = get_snps_from_segs(orig_Q(conflict_segs,:));
            conflict_genes = find_genes_from_snps(gene_gistic_rg(genes_in_chr),conflict_snps);
            new_sc(conflict_genes) = update_gene_gistic_scores(Z.dat{k},...
                            gene_gistic_rg(genes_in_chr(conflict_genes)),...
                            gene_gistic_ds,gene_gistic_res);
          end
        else
          % no conflict or no arbitration
          new_sc = sc;
        end
        
        if max(new_sc) >= score_thresh(k)
          if ~isempty(conflict_segs) && do_arbitration
            seg_assignment{k}(conflict_segs,:) = dist_weights;
          end
          if gene_gistic
            [mx gene_mi] = max(new_sc);
            top_genes = find(new_sc == mx);
            mi = gene_gistic_rg(genes_in_chr(gene_mi)).midsnp-chr_zero;
            rg = [];
            for i=1:length(top_genes)
              rg = [rg; gene_gistic_rg(genes_in_chr(top_genes(i))).snps-chr_zero];
            end
            rg = unique(rg);
          else
            [mx mi] = max(new_sc);
            if isempty(intersect(mi,rg))
              warning('Arbitrated peel-off has shifted peak!');
            end
            rg = find(new_sc == mx);
            [rg,mi] = fix_rg(rg,mi);
          end
          regs_on_chr = [regs_on_chr n];
        else
          % no longer any significant scores
          % shrink assignment table and exit
          seg_assignment{k} = seg_assignment{k}(:,1:n-1);
          n = n-1;
          loop = 0;
          break;
        end
 
        %% current peak is significant: create a regs entry for it and
        %% save some basic information about it

        regs{k}(n).chrn = ch; % sets chr number of current reg
        if gene_gistic
          [~,longest_top] = max(cellfun(@length,{gene_gistic_rg(genes_in_chr(top_genes)).snps}));
          qq = gene_gistic_rg(genes_in_chr(top_genes(longest_top))).snps-chr_zero;
          [~,jj] = find(y(qq,:));
          samples = unique(jj);
        else          
          samples = find(y(mi,:));
        end
        
        v=zeros(1,ns); 
        v(samples)=1; % sets samples with alteration to 1
        regs{k}(n).samples=v; 
        regs{k}(n).peak_st = min(rg)+chr_zero;
        regs{k}(n).peak_en = max(rg)+chr_zero;
        regs{k}(n).peak=round(0.5*(regs{k}(n).peak_st+regs{k}(n).peak_en));
        if gene_gistic
          % gene gistic q-value is from gene 
          regs{k}(n).qv = gene_gistic_rg(genes_in_chr(gene_mi)).q;
        else
          % non-gene-gistic q-value is from SNP
          regs{k}(n).qv = q{k}(chr_zero+mi); 
        end
        nsnps = length(Z.pos);
        if gene_gistic
          regs{k}(n).st = min(nsnps,max([0;find(orig_ads(1:mi)<score_thresh(k))])+chr_zero+1);
          right_half = (mi:length(orig_ads))+1;
          regs{k}(n).en = min([ right_half(find(orig_ads(mi:end) < ...
                                              score_thresh(k)))+chr_zero+1 chr_max]);

          regs{k}(n).gene_symbs = unique({gene_gistic_rg(genes_in_chr(top_genes)).symb});
          regs{k}(n).orig_score = orig_sc(gene_mi);
          regs{k}(n).iter_score = mx;
        else
          % find first snp on chr scoring above score_thresh
          regs{k}(n).st = min(nsnps,max([0; find(orig_sc(1:mi)<score_thresh(k))])+chr_zero+1);
          % find last snp on chr scoring above thresh  
          right_half = (mi:length(sc))+1;
          regs{k}(n).en = min([ right_half(find(orig_sc(mi:end) < ...
                                              score_thresh(k)))+chr_zero+1 chr_max]);
          regs{k}(n).orig_score = orig_sc(mi);
          regs{k}(n).iter_score = mx;
        end
        
        % Recalculate score and find new maximum
        
        Z.dat{k} = Znew;
        Q = Qnew;
        x = Z.dat{k}(in_chr,:);
        y=(x>0);
        
        % update score and find new maximum peak 
        if gene_gistic
          peel_off_snps = get_snps_from_segs(Qrm);
          peel_off_genes = find_genes_from_snps(gene_gistic_rg(genes_in_chr),peel_off_snps);
          sc(peel_off_genes) = update_gene_gistic_scores(Z.dat{k},...
                                    gene_gistic_rg(genes_in_chr(peel_off_genes)),...
                                    gene_gistic_ds, gene_gistic_res);
          [mx gene_mi] = max(sc);
          top_genes = find(sc == mx);
          mi = gene_gistic_rg(genes_in_chr(gene_mi)).midsnp-chr_zero;
          rg = [];
          for i=1:length(top_genes)
            rg = [rg; gene_gistic_rg(genes_in_chr(top_genes(i))).snps-chr_zero];
          end
          rg = unique(rg);
        else
          sc = nanmean(x,2);
          [mx mi] = max(sc); % find maximum score and max snp
          rg=find(sc==mx); % find entire peak region where sc = mx
   
          % pick larger of any noncontiguous regions
          [rg,mi] = fix_rg(rg,mi);
          
        end
        
        
      end % end of main peak finding loop
      
      %% there are no more significant peaks on this chromosome
                     
      % Now, go through and calculate post-arbitration scores, resid-qvs,
      % and robust regions for each peak 
            
      % First, find unassigned segments.   These are the 'passenger' segs
      % and should be included in the score for all peaks.
      
      chr_segs = find(orig_Q(:,1) == ch);
      assigned_segs = find(sum(seg_assignment{k}(:,regs_on_chr),2));
      unassigned_segs = setdiff(chr_segs,assigned_segs);
      
      verbose('Computing background score...',20)
      
      if use_segarray
        %! use SegArray.addSegments to add unassigned segments to empty background
        bg_Z = SegArray.zeros(nnz(in_chr),ns);
        bg_Z = addSegments(bg_Z,orig_Q(unassigned_segs,2)-chr_zero, ...
                              orig_Q(unassigned_segs,3)-chr_zero, ...
                              orig_Q(unassigned_segs,5), ...
                              orig_Q(unassigned_segs,4) );
      else % ~use_segarray
        % add unassigned segments one at a time to initially empty background
        bg_Z = zeros(nnz(in_chr),ns);
        for j=1:length(unassigned_segs)
          cur_seg = orig_Q(unassigned_segs(j),:);
          bg_Z(cur_seg(2)-chr_zero:cur_seg(3)-chr_zero,cur_seg(5)) = ...
              bg_Z(cur_seg(2)-chr_zero:cur_seg(3)-chr_zero,cur_seg(5)) + ...
              cur_seg(4);
        end
      end

      % loop across regions of the chromosome
      for i=1:length(regs_on_chr)
        verbose('Computing score for peak %d of %d on chromosome',30,[i length(regs_on_chr)]);
        cur_reg = regs_on_chr(i);
        add_segs = find(seg_assignment{k}(:,cur_reg));
        if use_segarray
          % use SegArray.addSegments to compute score
          cur_Z = addSegments(bg_Z,orig_Q(add_segs,2)-chr_zero, ...
                              orig_Q(add_segs,3)-chr_zero, ...
                              orig_Q(add_segs,5), ...
                              seg_assignment{k}(add_segs,cur_reg) ); %!
        else % not using segarray
          cur_Z = bg_Z;
          for j=1:length(add_segs)
            cur_seg = orig_Q(add_segs(j),:);
            amt_to_add = seg_assignment{k}(add_segs(j),cur_reg);
            cur_Z(cur_seg(2)-chr_zero:cur_seg(3)-chr_zero,cur_seg(5)) = ...
                cur_Z(cur_seg(2)-chr_zero:cur_seg(3)-chr_zero,cur_seg(5)) + ...
                amt_to_add;
          end
        end
        
        if ~gene_gistic
          cur_ads = nanmean(cur_Z,2);
          %regs{k}(cur_reg).score = sum(seg_assignment{k}(:,cur_reg))/ ...
          %  size(Z.dat{k},2);
          regs{k}(cur_reg).score = cur_ads(regs{k}(cur_reg).peak- ...
                                           chr_zero);
          regs{k}(cur_reg).ads = cur_ads;
        else
          cur_sc = update_gene_gistic_scores(cur_Z, ...
                                             gene_gistic_rg(genes_in_chr), ...
                                             gene_gistic_ds, ...
                                             gene_gistic_res,chr_zero);
          peak_genes = [];

          % accumulate indices to peak genes
          for j=1:length(regs{k}(cur_reg).gene_symbs)
            peak_genes = unique([peak_genes;
                                strmatch(regs{k}(cur_reg).gene_symbs{j},{gene_gistic_rg(genes_in_chr).symb})]);
          end
          regs{k}(cur_reg).score = max(cur_sc(peak_genes));
          regs{k}(cur_reg).ads = cur_sc;
          regs{k}(cur_reg).peak_genes = peak_genes;
        end
        
        regs{k}(cur_reg).resid_qv = exp(interp_pwl(sorted_ads,sorted_q,regs{k}(cur_reg).score));
        
        %% Calculate robust region start and end 
        cur_mx = regs{k}(cur_reg).score;
        if ~gene_gistic
          left_half = 1:regs{k}(cur_reg).peak_st-chr_zero;
          right_half = regs{k}(cur_reg).peak_en-chr_zero:nnz(in_chr);
          regs{k}(cur_reg).robust_reg_st = max([1;find(cur_ads(left_half)<(cur_mx- score_thresh(k)))])+chr_zero; 
          regs{k}(cur_reg).robust_reg_en = min([ right_half(find(cur_ads(right_half)< (cur_mx-score_thresh(k))))+chr_zero chr_max]); 
        else
          left_half = 1:min(regs{k}(cur_reg).peak_genes);
          right_half = max(regs{k}(cur_reg).peak_genes):length(genes_in_chr);
          regs{k}(cur_reg).robust_reg_st_gene = max([0 find(cur_sc(left_half)< (cur_mx - score_thresh(k)))]);
          regs{k}(cur_reg).robust_reg_en_gene = ...
              min([right_half(find(cur_sc(right_half)< (cur_mx-score_thresh(k)))) Inf]); 
        end
                  
        if k==1
          disp(['Amplification peak at ' genomic_location(Z,{regs{k}(cur_reg).peak_st:regs{k}(cur_reg).peak_en})]);
        else
          if ~gene_gistic
            disp(['Deletion peak at ' genomic_location(Z, ...
                                                       {regs{k}(cur_reg).peak_st:regs{k}(cur_reg).peak_en})]);
          else
            disp(['Deletion peak at genes: ' ...
                  regs{k}(cur_reg).gene_symbs{1}]);
          end
        end
        disp(['q-value = ' num2str(regs{k}(cur_reg).qv) ' & residual q-value = ' num2str(regs{k}(cur_reg).resid_qv)]);
        
      end % matches for i=1:length(regs_on_chr)
      
    end %% end of ch loop
    
  end %% end of k=1:2 loop
  
    
    
          
    
  
                    
  
