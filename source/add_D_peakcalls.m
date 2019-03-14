function [D,supidx]=add_D_peakcalls(D,calls,prefix,regs,cyto,pvs)


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

if ~exist('prefix','var') || isempty(prefix)
  prefix={'',''};
end
if ischar(prefix)
  prefix={prefix,prefix};
end

if exist('cyto','var')
  D=add_cyto(D,cyto);
end

if ~exist('pvs','var') || length(pvs)==1
  for k=1:length(regs)
    pvs{k}=zeros(size(D.dat,1),1);
  end
end

acc={'AP','DP'};
desc={'Amplification Peak #','Deletion Peak #'};
supidx = cell(1,length(regs));
for k=1:length(regs)
  if isempty(regs{k})
    continue;
  end
  ncalls = size(calls{k},1);
  barsep = repmat('| ',ncalls,1);
  peakacc=[ repmat([prefix{k} acc{k}],size(calls{k},1),1) num2str((1:size(calls{k},1))','%d')];
  if ~exist('regs','var')
    peakdesc=[ repmat([prefix{k} desc{k}],size(calls{k},1),1) ...
               num2str((1:size(calls{k},1))','%d')];
  else
    % make wide peak column
    if isfield(regs{k}(1),'peak_wide_st')
      wide_peak_str=[ repmat(['| wide_peak '],size(calls{k},1),1) ...
                     strvcat(genomic_location(D,mat2cell([cat(1,regs{k}.peak_wide_st) cat(1,regs{k}.peak_wide_en)],...
                                                        ones(length(regs{k}),1),2),cyto,1)) ];
    else
      b_wide_peak_str=[];
    end
    
    % make bootstrap wide peak column (currently not used)
    if isfield(regs{k}(1),'b_wide_peak_st')
      b_wide_peak_str=[ repmat(['| bootstrp_wide_peak '],size(calls{k},1),1) ...
                     strvcat(genomic_location(D,mat2cell([cat(1,regs{k}.b_wide_peak_st) cat(1,regs{k}.b_wide_peak_en)],...
                                                        ones(length(regs{k}),1),2),cyto,1)) ];
    else
      b_wide_peak_str=[];
    end
    
    % make residual q-value column
    if isfield(regs{k}(1),'resid_qv')
      resid_qv_str=[ repmat(['| resid q-val '],size(calls{k},1),1) ...
                      strvcat(num2str(cat(1,regs{k}.resid_qv))) ];
    else
      resid_qv_str=[];
    end
    
    % make q-value column
    if isfield(regs{k}(1),'qv')
      orig_qv_str = [barsep num2str([regs{k}.qv]')];
    else
      orig_qv_str = [barsep num2str(pvs{k}(cat(1,regs{k}.peak_st)))]; %! untested
    end
    
    if isfield(regs{k}(1),'focal')
        focal_broad_options={'none','focal','broad','both'};
        tmp_f=cat(1,regs{k}.focal);
        tmp_b=cat(1,regs{k}.broad);
        focal_broad_st=strvcat(focal_broad_options(tmp_f+tmp_b*2+1));
    elseif isfield(regs{k}(1),'broad')
        focal_broad_options={'none','broad'};

        tmp_b=cat(1,regs{k}.broad);
        focal_broad_st=strvcat(focal_broad_options(tmp_b+1));
    else
        focal_broad_st=[];
    end

    peakdesc=[ peakacc ... 
               [repmat(['| region '],size(calls{k},1),1) ...
                strvcat(genomic_location(D,mat2cell([cat(1,regs{k}.st) cat(1,regs{k}.en)],...
                                                    ones(length(regs{k}),1),2),cyto,1)) ...
                repmat('| peak ',size(calls{k},1),1) ...
                strvcat(genomic_location(D,mat2cell([cat(1,regs{k}.peak_st) cat(1,regs{k}.peak_en)],...
                                                    ones(length(regs{k}),1),2),cyto,1)) ...
                wide_peak_str ...
                barsep b_wide_peak_str ...
                barsep strvcat(cyto(D.cyton(cat(1,regs{k}.peak))).name) ...
                orig_qv_str ...
                resid_qv_str ...
                barsep focal_broad_st ...
               ] ];
  end
  [D,supidx{k}]=add_D_sup(D,peakacc,peakdesc,calls{k},'cols');
end




