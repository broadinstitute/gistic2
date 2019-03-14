function gene_calls(run_dir,params,refgene)
%GENE_CALLS create copy number alteration calls for genes on all samples
%
%     gene_calls(RUN_DIR,PARAMS,RG_FILE)
%     gene_calls(RUN_DIR,EXT,RG_FILE)  [older form]
%
% This is a file-oriented function which assigns copy number alteration 
% values to genes for every sample using the file 'all_data_by_genes.txt'
% as input. The assigned values are:
%
%     -2 => high level deletion
%     -1 => significant deletion
%      0 => no significant copy number alteration
%      1 => significant amplification
%      2 => high level amplification
%
% The thresholds for high level amplifications are determined
% sample-by-sample using sample cutoffs generated either from the
% 'broad_values_by_arm.txt' output file (if available) or
% 'broad_data_by_genes.txt' (if arm values are not available).
% The data are output to 'all_thresholded.by_genes.txt' and the sample
% cutoffs used are output to 'sample_cutoffs.txt'.
%
% The RUN_DIR parameter names the working directory where the output files
% (some of which are inputs to this function) are located. The PARAMS
% structure contains information about thresholds and file naming
% conventions:
%       PARAMS.t_amp and .t_del are the amplification and deletion
%     thresholds used in the GISTIC analysis
%       PARAMS.fname and PARAMS.ext name the prefix and extension used in
%       the file names
% The second argument for older form of GENE_CALLS was just the extension
% (no file name will be used). RG_FILE is either a path to the file
% that genome features, or the contents of that file loaded into a struct

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


  if isstruct(params)
    % if second argument is params structure
    ext = params.ext;
    fname = params.fname;
    t_amp = params.t_amp;
    t_del = params.t_del;
    genepattern = params.genepattern;
  else
    % compatible behavior - second argument was extension
    ext = params;
    load([run_dir 'ts' ext '.mat']); % read t_amp, t_del from file
    fname = '';
    genepattern = false;
  end

  pathbase = [run_dir fname]; % path + input file base name

  %% Step 1: Create sample_cutoffs
  verbose('Calculating sample cutoffs...',30);
%! load(cyto_file)
  if exist([pathbase 'broad_values_by_arm' ext '.txt'],'file')
    B=read_by_arms_file([pathbase 'broad_values_by_arm' ext '.txt']);
  else
    B=read_by_gene_file([pathbase 'broad_data_by_genes' ext '.txt'], ...
                        refgene);
  end
  arm_min=min(B.dat,[],1);
  arm_max=max(B.dat,[],1);
  cutoff_min=arm_min-t_del; % from dynamic cutoffs
                            % cutoff_max=max(arm_max,-arm_min)+t_amp; % from dynamic cutoffs
  cutoff_max=arm_max+t_amp; % going back to original method
  
  % Write sample cutoff file
  verbose('Writing sample cutoff file...',30);
  f=fopen([pathbase 'sample_cutoffs' ext '.txt'],'w');
  fprintf(f,'# Using amp_thresh=%4.2f, del_thresh=%4.2f\n',t_amp,t_del);
  fprintf(f,'Sample\tLow\tHigh\n');
  for i=1:length(arm_max)
    fprintf(f,'%s\t%f\t%f\n',B.sdesc{i},cutoff_min(i),cutoff_max(i));
  end
  fclose(f);
    
  %% Create thresholded data
  verbose('Thresholding data...',30);
  tbl=read_table([pathbase 'sample_cutoffs' ext '.txt'],'%s%f%f',char(9),2,'commentstyle','shell');
  A=read_by_gene_file([pathbase 'all_data_by_genes' ext '.txt'],refgene);                      
  [~,m1,m2]=match_string_sets_hash(tbl.dat{1},A.sdesc);
  tmp=cat(2,tbl.dat{2:3})';
  A.ssupdat=nan(2,size(A.dat,2));
  A.ssupdat(1:2,m2)=tmp(:,m1);
  A.ssupacc={'Low','High'};
  A.ssupdesc={'Low','High'};
%!unused  OA=order_by_pos(A);
  T=A;
  T.dat=zeros(size(A.dat));
  for i=1:size(T.dat,2)
    T.dat(A.dat(:,i)<A.ssupdat(1,i),i)=-2;
    T.dat(A.dat(:,i)>=A.ssupdat(1,i) & A.dat(:,i)<-t_del,i)=-1;
    T.dat(A.dat(:,i)>t_amp & A.dat(:,i)<=A.ssupdat(2,i),i)=1;
    T.dat(A.dat(:,i)>A.ssupdat(2,i),i)=2;
  end
%!  OT=order_by_pos(T);

  %% write all_thresholded.by_genes.txt
  verbose('Writing thresholded gene calls...',30);
  % matlab output
  if ~genepattern
    save([pathbase 'all_thresholded.by_genes' ext '.mat'],'T');
  end
  % tab-delimited text output
  f=fopen([pathbase 'all_thresholded.by_genes' ext '.txt'],'w');
  fprintf(f,'Gene Symbol\tLocus ID\tCytoband');
  for j=1:size(T.dat,2)
    fprintf(f,'\t%s',T.sdesc{j});
  end
  fprintf(f,'\n');
  for i=1:size(T.dat,1);
    modi(i,100)
    fprintf(f,'%s\t%s\t%s',T.gacc{i},T.gdesc{i},T.cyto{i});
    for j=1:size(T.dat,2)
      fprintf(f,'\t%d',T.dat(i,j));
    end
    fprintf(f,'\n');
  end
  fclose(f);
