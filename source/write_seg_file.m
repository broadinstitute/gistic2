function write_seg_file(fname,D,convert_to_raw,report_nan_segments)
%WRITE_SEG_FILE write copy number data out as a segmented data file

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


% defaults for optional arguments
if ~exist('convert_to_raw','var') || isempty(convert_to_raw)
  convert_to_raw = false;
end
if ~exist('report_nan_segments','var') || isempty(report_nan_segments)
  report_nan_segments = false;
end

% data may or may not be in log space
if isfield(D,'islog')
  islog = D.islog;
else
  % assume islog if unknown
  islog = true;
end

% create segmented data
rl=runlength(D.dat,D.chrn);

if ~iscell(rl)
    rl = {rl};
end

if iscell(D.sdesc)
    D.sdesc = char(D.sdesc);
end

f=fopen(fname,'w');
% print appropriate header
if convert_to_raw | ~islog;
  fprintf(f,'%s\t%s\t%s\t%s\t%s\t%s\n','Sample','Chromosome','Start.bp','End.bp','Num.Markers','Seg.CN');
else
  fprintf(f,'%s\t%s\t%s\t%s\t%s\t%s\n','Sample','Chromosome','Start.bp','End.bp','Num.Markers','Log2.Ratio');
end
% output data sample-by-sample
for i=1:length(rl)
  rli = rl{i};
  if report_nan_segments
    idx=(1:size(rli,1))';
  else
    idx=find(~isnan(rli(:,3)));
  end
  % segment-by-segment within sample
  for j=idx'
    if convert_to_raw
      fprintf(f,'%s\t%s\t%d\t%d\t%d\t%f\n',char(D.sdesc(i,:)),num2chromosome(D.chrn(rli(j,1))),...
              D.pos(rli(j,1)),D.pos(rli(j,2)),rli(j,2)-rli(j,1)+1,2^(rli(j,3)+1));
    else
      fprintf(f,'%s\t%s\t%d\t%d\t%d\t%f\n',char(D.sdesc(i,:)),num2chromosome(D.chrn(rli(j,1))),...
              D.pos(rli(j,1)),D.pos(rli(j,2)),rli(j,2)-rli(j,1)+1,rli(j,3));
    end
  end
end
fclose(f);
