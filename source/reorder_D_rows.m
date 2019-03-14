function D=reorder_D_rows(D,varargin)
%REORDER_D_ROWS  Reduce data set to selected rows.
%
%   D = REORDER_D_ROWS(D,INDICES) unpacks fields from data structure D and
%   repacks using a subset of the data given by the column indices INDICES.
%
%   D = REORDER_D_ROWS(D,GSUPACC_NAME,VAL) filters the data structure for
%   rows whos GSUPACC_NAME in .gsupacc matches VAL in .gsupdat.
%
%   D = REORDER_D_ROWS(D,IDX,'allmem') Same as reorder_D_ROWS(d,idx).
%   %   (Case supported for datastruct compatibility).
%
%   D = REORDER_D_ROWS(D,SUPACC,VAL,'allmem') Same as
%   reorder_D_ROWS(d,supacc,val) (Case supported for datastruct
%   compatibility).

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


%   In both usages, data structure history is updated.
%
%   Example:  D = reorder_D_rows(D,[1:10:size(D.dat,1)]) decimates the data 
%                 structure, returning a new D with every 10th row (SNP) 
%                 included.
%
%             D = reorder_D_rows(D,'Interesting',1) returns a new data
%             structure consisting only of SNPs with .gsupacc  of 'Interesting'
%             set to 1 in the .gsupdat field.
%
%   History
%       -- 09 Oct 07  --- Added catch for if field exists but is empty at
%       line 123.
%
%`      -- 17 Oct 07  --- Added .plat to reorder fields.
%---
% $Id$
% $Date$
% $LastChangedBy$
% $Rev$

if length(varargin)==0
  error('missing argument');
elseif length(varargin)==1
  idx=varargin{1};
else
  idx=find_D_rows(D,varargin{1},varargin{2});
end

%% Parse varargin

if isempty(varargin)
  error('missing argument')
elseif length(varargin)>=2 && ischar(varargin{1})
  idx=find_D_rows(D,varargin{1},varargin{2});
else
  idx = varargin{1};
end

donefields = {};

%% Get origsize1

if isfield(D,'dat')
    D=add_history(D,mfilename,varargin{:});
    origsize1 = size(D.dat,1);
 
elseif isfield(D,'chrn')
    origsize1 = length(D.chrn);
elseif isfield(D,'pos')
    origsize1 = length(D.pos);
elseif isfield(D,'marker')
    origsize1 = length(D.marker);
elseif isfield(D,'affy_calls')
    origsize1 = size(D.affy_calls,1);
end

%Update gdesc field


 %Update gdesc field

if isfield(D,'gdesc')  
  if iscell(D.gdesc)
    D.gdesc = D.gdesc(idx);
  else
    D.gdesc = deblank(D.gdesc(idx,:));
  end
  donefields = union(donefields,{'gdesc'});
end

if isfield(D,'gsymb')  %Update gsymb field
  if iscell(D.gsymb)
    D.gsymb = D.gsymb(idx);
  else
    D.gsymb = deblank(D.gsymb(idx,:));
  end
  donefields = union(donefields,{'gsymb'});
end

if isfield(D,'gacc')   %Update gacc field
  if iscell(D.gacc)
    D.gacc = D.gacc(idx);
  else
    D.gacc = deblank(D.gacc(idx,:));
  end
  donefields = union(donefields,{'gacc'});
end

if isfield(D,'id')      %Update id field
  if iscell(D.id)
    D.id = D.id(idx);
  else
    D.id = deblank(D.id(idx,:));
  end
  donefields = union(donefields,{'id'});
end

if isfield(D,'symb')    %Update symb field
  if iscell(D.symb)
    D.symb = D.symb(idx);
  else
    D.symb = deblank(D.symb(idx,:));
  end
  donefields = union(donefields,{'symb'});
end

if isfield(D,'affy_call')   %Update affy_call field
  D.affy_call = D.affy_call(idx,:);
  donefields = union(donefields,{'affy_call'});
end

if isfield(D,'adat')        %Update adat field
  D.adat = D.adat(idx,:,:);
  donefields = union(donefields,{'adat'});
end

if isfield(D,'gsupdat')     %Update gsupdat field
  if size(D.gsupdat,1)==size(D.dat,1) % rows match
    warning('recommend using a column per row in dat'); %#ok<WNTAG>
    D.gsupdat=D.gsupdat(idx,:);
  else
    D.gsupdat=D.gsupdat(:,idx);    
  end
  donefields = union(donefields,{'gsupdat'});
end

% Update remaining vector fields
vector_fields={'cyto','gorigidx','marker','chr','chrn','cM','pos','score','grg','start','end','ll','cyton','cyto_stain', ...
               'armn','chrloc','refgene','refgene_idx','grange','chrarmn','n_collapse','plat'};
if isfield(D,'gene_fields')
  vector_fields = unique([vector_fields D.gene_fields]);
end
for j=1:length(vector_fields)
  field = vector_fields{j};
  if isfield(D,field) && ~ismember(field,donefields)
    x = D.(field);
    if isvector(x)
      x = x(idx);
    else
      x = x(idx,:);
    end
    D.(field) = x;
    donefields = union(donefields,{field});
  end
end

% Update remaining matrix fields
matrix_fields={'dat','smooth','sm1','sm2','sm3','sm2j','raw','cbs','cbs_fixed','flag','hmm','affy_calls','ref','ads','qv','fxa','fxa_extra','pvs'};
if isfield(D,'matrix_fields')
    matrix_fields = unique([matrix_fields D.matrix_fields]);
end
for j=1:length(matrix_fields)
    field = matrix_fields{j};
    if isfield(D,field)
        x = D.(field);
        if ~isempty(x)
            x = x(idx,:);
            D.(field) = x;
        end
        donefields = union(donefields,{field});
    end
end

%Now take care of fields that aren't hard-coded here

if exist('origsize1','var')
    fields = fieldnames(D);
    remfields = setdiff(fields,[matrix_fields vector_fields {'dat','gsesc','gsymb','gacc','id','symb','affy_call','adat','gsupdat'}]);

    for j = 1:size(remfields)
        fld = remfields{j};
        if size(D.(fld),1) == origsize1
            D.(fld) = D.(fld)(idx,:);
        else
            continue
        end
    end
end
