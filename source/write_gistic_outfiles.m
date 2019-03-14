function write_gistic_outfiles(D,base_dir,ext,cyto,regs,ts,pvs,p,q,ads,...
                               rg,add_vals,score_thresh,find_broad,...
                               genepattern,partial_hits,islog,fname)
%WRITE_GISTIC_OUTFILES Write tab-delimited GISTIC output files

% GISTIC software version 2.0
% Copyright (c) 2011, 2016 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


  if ~exist('find_broad','var') || isempty(find_broad)
    find_broad = 1;
  end

  if ~exist('genepattern','var')
    genepattern = [];
  end

  if ~exist('partial_hits','var')
    partial_hits = [];
  end
  
  if ~exist('islog','var') || isempty(islog)
    islog = 1;
  end
  
  if ~exist('fname','var')
    fname = '';
  end
    
  all_lesion_fname = [base_dir fname 'all_lesions' ext '.txt'];
  disp(['Writing all_lesions file to: ' all_lesion_fname]);
  start_at=[];
  no_call_thresh=ts;
  
  min_cutoff = zeros(1,2);
  for k=1:2
    min_val = min(ads{k}(q{k} < 1));
    if ~isempty(min_val)
      min_cutoff(k)=min_val;
    else
      min_cutoff(k) = 0;
    end
  end
  if ~exist('p_arm','var')
    p_arm=struct('method','scorecutoff','ads',{ads},'p_arm',0.5,'score_thresh',score_thresh,...
                 'score_thresh_focal',score_thresh,'min_cutoff',min_cutoff);
  end		
 
  % focal_broad output table (Note: 'find_broad' has been hardwired to 0
  % from the caller, run_focal_gistic for some time)
  if find_broad
    regs1=find_broad_regs(D,cyto,regs,p_arm);
    f=fopen([base_dir fname 'focal_broad' ext '.txt'],'w');
    ampdel={'Amp','Del'};
    for k=1:2
      for i=1:length(regs1{k})
        [st,~,~]=genomic_location(D,{[regs1{k}(i).st regs1{k}(i).en]},cyto,1,0);
        fprintf(f,'%s\n',[ ampdel{k} num2str(i) ')' st ' focal:' num2str(regs1{k}(i).focal) ' broad:' num2str(regs1{k}(i).broad)]);
      end
    end
    fclose(f);
  end
  
  % output the all_lesions table
  make_all_lesions_file(all_lesion_fname,D,regs,pvs,cyto,start_at,ts,no_call_thresh,add_vals,p_arm,find_broad,islog);

  %% output gene-annotated peaks {amp|del}_genes and deprecated table_{amp|del}
  
  % add .symbol field to reference genome
  % (NOTE: .symbol was supposed to add marks to .symb matching genes specified in 
  %  an annotation file, but this has been hardwored off in code below
  disp('Adding gene annotations'); 
  rg=add_chrn(rg);
  annot_file=[];
  mark=[];
  rg=add_annotation(rg,annot_file,mark);
  calls=call_regs(D,regs,{ts});
  
  % call function to write the tables
  disp('Writing output gene table files'); 
  A=[];
  u95ll=[];
  genetables(D,A,u95ll,rg,cyto,regs,calls,ts,ext,[],base_dir,genepattern,partial_hits,fname);
  fclose('all');
