function [CL21x,supids]=make_all_lesions_file(fname,CL21,regs,pvs,cyto,start_at,ts,no_call_thresh,add_vals,p_arm,find_broad,islog)
% [CL21x,supids]=make_all_lesions_file(fname,CL21,regs,pvs,cyto,start_at,ts,no_call_thresh,add_vals,p_arm)
%
%  TO MODIFY ALL LESIONS FILE CHANGE SUBFUNCTION 'FORMATFILE' BELOW

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


%% PUT DATA TO WRITE IN CL21X.SUPDAT

% add START line

tsx2 = [0.9 1.3];

if ~exist('find_broad','var') || isempty(find_broad)
  find_broad = 1;
end

if ~exist('ts','var')
  ts=[0.3 0.3];
end

if ~exist('islog','var') || isempty(islog)
  islog = 1;
end

if ~exist('no_call_thresh','var')
  no_call_thresh=ts;
end

if ~exist('add_vals','var')
  add_vals=0;
end

if ~exist('p_arm','var')
  p_arm=0.5;
end

if isempty(regs)
  CL21x=CL21;
  supidx=[];
  return;
end

v=zeros(1,size(CL21.dat,2));
CL21=add_D_sup(CL21,'START','START',v,'cols');

if exist('start_at','var') && ~isempty(start_at) && size(CL21.supdat,1)>1
  idx=1:size(CL21.supdat,1);
  idx=[ idx(1:(length(idx)+start_at-1)) length(idx) idx((length(idx)+start_at):(end-1))];
  CL21=reorder_D_sup(CL21,'cols',idx);
end

% generate broad regions

if find_broad
  regs=find_broad_regs(CL21,cyto,regs,p_arm);
end

for k = 1:length(regs)
    if ~isempty(regs{k})
      peakstarts = [regs{k}.peak_st];
      peakchrn = [regs{k}.chrn];
      posn = peakstarts + 1e11.*peakchrn;
      [dum,i] = sort(posn);
      regs{k} = regs{k}(i);
    end
end

    
   
if find_broad
  [CL21x,supids{1}]=add_polysomy_by_reg(CL21,regs,ts,0.5,no_call_thresh);
else
  CL21x = CL21;
end

% add x1 x2 any calls
[calls,vals]=call_regs(CL21,regs,{ts},{no_call_thresh});

calls_x2 = {[],[]};
calls_x1 = {[],[]};

for k=1:length(regs)
  if isempty(regs{k})
    continue
  end
  switch k
   case 1
    calls_x2A=call_regs(CL21,regs,{tsx2});
    calls_x2{1}=calls_x2A{1};
    calls_x1{1}=calls{1}-calls_x2{1};
   case 2
    calls_x2D=call_regs(CL21,regs,{tsx2});
    calls_x2{2}=calls_x2D{2};
    calls_x1{2}=calls{2}-calls_x2{2};
  end
end
[CL21x,supids{2}]=add_D_peakcalls(CL21x,calls,'Any-',regs,cyto,pvs);
[CL21x,supids{3}]=add_D_peakcalls(CL21x,calls_x2,'x2-',regs,cyto,pvs);
[CL21x,supids{4}]=add_D_peakcalls(CL21x,calls_x1,'x1-',regs,cyto,pvs);
if add_vals
  [CL21x,supids{4}]=add_D_peakcalls(CL21x,vals,'log2rat-',regs,cyto,pvs);
end



%% FORMAT DATA
[towrite,formatstring] = formatfile(CL21x,ts,tsx2,islog);

%% WRITE FILE
fid = fopen(fname,'w');
fprintf(fid,formatstring,towrite{:});
fclose(fid);

end

%%  
function [towrite,formatstring] = formatfile(CL21x,ts,tsx2,islog)

if ~iscell(CL21x.sdesc)
    CL21x.sdesc = cellstr(CL21x.sdesc);
end



%%%%%%%MODIFY HERE TO CHANGE COLUMN NUMBERS 
cols.un = 1;  %unique name
cols.d = 2;  %descriptor
cols.wl = 3;  %wide peak limit
cols.pl = 4;  %peak limit
cols.rl = 5;  %region limit
cols.qv = 6;  %q vaules
cols.rq = 7; %residual q value
cols.bf = 8;  %broad/focal
cols.at = 9;  %amplitude threshold
cols.dat = 10:(10+length(CL21x.sdesc)-1); %sample data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%MODIFY HERE TO CHANGE COLUMN HEADINGS
header{cols.un} = 'Unique Name';
header{cols.d} = 'Descriptor';
header{cols.wl} = 'Wide Peak Limits';
header{cols.pl} = 'Peak Limits';
header{cols.rl} = 'Region Limits';
header{cols.qv} =  'q values';
header{cols.rq} = 'Residual q values after removing segments shared with higher peaks';
header{cols.bf} = 'Broad or Focal';
header{cols.at} = 'Amplitude Threshold';
header(cols.dat) = CL21x.sdesc';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% FORMAT PEAK DATA with 0,1,2 values  (peakdata)

anysupaccidx = strmatch('Any-',CL21x.supacc);  %Get supdat indices for peakdata

%TEXTFROMSUP returns structure with fields same as COLS and containing text corresponding to
%fields of cols
peakstext = textfromsup(CL21x,anysupaccidx);  
peakstext.un = regexprep(cellstr(CL21x.supacc(anysupaccidx,:)),'Any-','');

%%%%%%%%% INSERT AMPLITUDES THRESHOLD STRING
ampthreshstr = ['0: t<' num2str(ts(1)) '; 1: ' num2str(ts(1)) '<t< ' num2str(tsx2(1)) '; 2: t>' num2str(tsx2(1))];
delthreshstr = ['0: t>-' num2str(ts(2)) '; 1: -' num2str(ts(2)) '>t> -' num2str(tsx2(2)) '; 2: t< -' num2str(tsx2(2))];
peakstext.at(strmatch('AP',peakstext.un)) = {ampthreshstr};
peakstext.at(strmatch('DP',peakstext.un)) = {delthreshstr}; 

%%%%%%%%%  FORMAT DATA FOR SAMPLES
arrayvals = CL21x.supdat((anysupaccidx + length(anysupaccidx)),:) + CL21x.supdat(anysupaccidx,:);
arrayvals = mat2cell(arrayvals,ones(length(anysupaccidx),1),ones(length(CL21x.sdesc),1));
arrayvals = reshape(arrayvals,1,numel(arrayvals));
arrayvals = cellfun(@num2str,arrayvals,'UniformOutput',0);
arrayvals = reshape(arrayvals,length(anysupaccidx),length(CL21x.sdesc));
peakstext.dat = arrayvals;

%%%%%%%%%  EXPAND PEAK NAMES
ampsname = 'Amplification Peak ';
delsname = 'Deletion Peak ';
peakstext.un = regexprep(peakstext.un,'AP',ampsname);
peakstext.un = regexprep(peakstext.un,'DP',delsname);

%%%%%%%%%  GENERATE CELL TO WRITE PEAK DATA
for fl = fieldnames(cols)'
    peakswrite(1:length(anysupaccidx),cols.(char(fl))) = peakstext.(char(fl)); %#ok<AGROW>
end


%% FORMAT LOG2 RATIOS


supaccidx = strmatch('log2rat',CL21x.supacc);
if ~isempty(supaccidx)
logstext = textfromsup(CL21x,supaccidx);

%%%%%%%%%  FORMAT DATA FOR SAMPLES
log2vals = CL21x.supdat(supaccidx,:) ;
log2vals = mat2cell(log2vals,ones(length(peakstext.un),1),ones(length(CL21x.sdesc),1));
log2vals = reshape(log2vals,1,numel(log2vals));
log2vals = cellfun(@num2str,log2vals,'UniformOutput',0);
log2vals = reshape(log2vals,length(peakstext.un),length(CL21x.sdesc));
logstext.dat = log2vals;
%%%%%%%
if islog
  logstext.at = repmat({'Actual Log2 Ratio Given'},length(supaccidx),1);
else
  logstext.at = repmat({'Actual Copy Change Given'},length(supaccidx),1);
end

logstext.un = strcat(peakstext.un,' - CN values');

%%%%%%%%%  GENERATE CELL TO WRITE LOG2 DATA
for fl = fieldnames(cols)'
    logswrite(1:length(supaccidx),cols.(char(fl))) = logstext.(char(fl)); %#ok<AGROW>
end

else
    logswrite = [];
end


%% FORMAT BROAD DATA

% THE SUP INDICES FOR BROAD DATA
broadidx = find(~cellfun(@isempty,regexp(cellstr(CL21x.supacc),'^[\d]')));

if ~isempty(broadidx)
  % THE PEAK INDICES CORRESPONDING TO BROADIDX
  idx = strmatch('b',peakstext.bf);
  
  isbroadidx = anysupaccidx(idx);
  
  %GET DATA FROM PEAK DATA THAT ARE BROADS
  broadstext = textfromsup(CL21x,isbroadidx);  
  % broadstext = textfromsup(CL21x,broadidx);
  
  
  %%%% CHANGE UNIQUE NAME FOR BROADS
  namebroads = strcat(cellstr(CL21x.supacc(broadidx,:)),'- Contributing Broad Events');
  broadsampname = ampsname;
  broadsdelname = delsname;
  namebroads = regexprep(namebroads,'(\d*)AB',[broadsampname ' $1']);
  namebroads = regexprep(namebroads,'(\d*)DB',[broadsdelname ' $1']);
  broadstext.un = namebroads;
  
  %%%% Change descriptor format
  broadsdescriptor = regexp(broadstext.d,'^[\dXY]+[\D]','match');
  broadstext.d = [broadsdescriptor{:}]';
  
  %%%%  Add suplemental information under peak limits and wide limits
  broadstext.pl = regexprep(strcat('Broad Event Corresponding to ', peakstext.un(idx)')','to','to ');
  broadstext.wl = repmat({'Amplitude values represent median within significant region.'},size(broadidx,1),1);
  
  
  
  broadstext.at(find(~cellfun(@isempty,regexp(broadstext.un,broadsampname)))) = {['0: t<' num2str(ts(1)) '; 1: t>' num2str(ts(1))]}; %#ok<FNDSB>
  broadstext.at(find(~cellfun(@isempty,regexp(broadstext.un,broadsdelname)))) = {['0: t>-' num2str(ts(2)) '; 1: t<-' num2str(ts(2))]};  %#ok<FNDSB>
  
  
  broaddat = CL21x.supdat(broadidx,:);
  broaddat = mat2cell(broaddat,ones(length(broadidx),1),ones(length(CL21x.sdesc),1));
  broaddat = reshape(broaddat,1,numel(broaddat));
  broaddat = cellfun(@num2str,broaddat,'UniformOutput',0);
  broadstext.dat = reshape(broaddat,length(broadidx),length(CL21x.sdesc));
  
  
  
  
  %%%%%%%%%  GENERATE CELL TO WRITE BROAD DATA
  for fl = fieldnames(cols)'
      if ~isempty(broadstext.(char(fl)))
        broadswrite(1:length(broadidx),cols.(char(fl))) = broadstext.(char(fl)); %#ok<AGROW>
      end
  end
else
  broadswrite={};
end

formatstring = [repmat('%s\t',1,length(header)) '\n'];

towrite = [header reshape(peakswrite',1,numel(peakswrite)) reshape(logswrite',1,numel(logswrite)) ...
           reshape(broadswrite',1,numel(broadswrite))];
end


function [ftext] = textfromsup(CL21x,supaccidx)

supdescs = strcat(strvcat(CL21x.supdesc(supaccidx,:)),'|'); %#ok<VCAT>

text = textscan(strvcat(supdescs)','%s%s%s%s%s%s%s%s%s','Delimiter','|');%#ok<VCAT> %,'MultipleDelimsAsOne',1); %#ok<VCAT>



descriptor = text{6};

widepeaklimits = text{4};
widepeaklimits = regexprep(widepeaklimits,'wide_peak ','');
widepeaklimits = regexprep(widepeaklimits,',chr.+[)]',')');
widepeaklimits = regexprep(widepeaklimits,'(','(probes ');

peaklimits = text{3};
peaklimits = regexprep(peaklimits,'peak ','');
peaklimits = regexprep(peaklimits,',chr.+[)]',')');
peaklimits = regexprep(peaklimits,'(','(probes ');

regionlimits = text{2};
regionlimits = regexprep(regionlimits,'region ','');
regionlimits = regexprep(regionlimits,',chr.+[)]',')');
regionlimits = regexprep(regionlimits,'(','(probes ');

qvalues = text{7};

qvaluespo = text{8};
qvaluespo = regexprep(qvaluespo,'resid q-val','');

amplitudethreshold = repmat({[]},length(supaccidx),1);

broadfocal = text{9};



ftext.un = repmat({[]},length(supaccidx),1);  %unique name
ftext.d = descriptor;  %descriptor
ftext.wl = widepeaklimits;  %wide peak limit
ftext.pl = peaklimits;  %peak limit
ftext.rl = regionlimits;  %region limit
ftext.qv = qvalues;  %q vaules
ftext.rq = qvaluespo; %residual q value
ftext.bf = broadfocal;  %broad/focal
ftext.at = amplitudethreshold;  %amplitude threshold
ftext.dat = repmat({[]},length(supaccidx),length(CL21x.sdesc)); %sample data


end







