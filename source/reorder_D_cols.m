function D=reorder_D_cols(D,varargin)
%REORDER_D_COLS - Reorder the column indices (samples) of D.
%
%   DNEW = REORDER_D_COLS(DOLD,IDX) reorders the columns in data structure
%   D according to vector IDX such that DNEW.SDESC = DOLD.SDESC(IDX).
%
%   DNEW = REORDER_D_COLS(DOLD,SUPACC,VAL) reorders the columns in data
%   structure D by finding the columns that match SUPPACC = VAL in the
%   SUPDAT field.
%
%   DNEW = REORDER_D_COLS(DOLD,IDX,'allmem') Same as reorder_D_cols(dold,idx).  %   (Case supported for datastruct compatibility).
%
%   DNEW = REORDER_D_COLS(DOLD,SUPACC,VAL,'allmem')  
%   Same as reorder_D_cols(dold,supacc,val)  
%   (Case supported for datastruct compatibility).

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


%---
% $Id$
% $Date%
% $LastChangedBy$
% $Rev$

D=add_history(D,mfilename,varargin{:});

%% Parse varargin
if isempty(varargin)
  error('Missing Argument')
elseif length(varargin) >= 2 && ischar(varargin{1})
  idx=find_D_cols(D,varargin{1},varargin{2});
else
    idx = varargin{1};
end


%% Get origsize2

if isfield(D,'dat')
    D=add_history(D,mfilename,varargin{:});
    origsize2 = size(D.dat,2);
elseif isfield(D,'chrn')
    origsize2 = length(D.chrn);
elseif isfield(D,'pos')
    origsize2 = length(D.pos);
elseif isfield(D,'marker')
    origsize2 = length(D.marker);
elseif isfield(D,'affy_calls')
    origsize2 = size(D.affy_calls,2);
end

%% Reorder fields

if isfield(D,'sdesc')
  if iscell(D.sdesc)
    D.sdesc=D.sdesc(idx);
  else
    D.sdesc=deblank(D.sdesc(idx,:));
  end
end

if isfield(D,'sscale') && ~isempty(D.sscale)
  D.sscale=D.sscale(idx,:);
end

if isfield(D,'adat')
  D.adat=D.adat(:,idx,:);
end

vector_fields={'sup','scans','scale2vec','residx','origidx','gcm_name','cbs_rl','sis','medians',...
               'peaks','joins','final_dist','sscaling'};
if isfield(D,'sample_fields')
    vector_fields = unique([vector_fields D.sample_fields]);
end
for j=1:length(vector_fields)
    if isfield(D,vector_fields{j})
        if ~isempty(D.(char(vector_fields{j})))
            x=D.(vector_fields{j});
            if isvector(x)
                x=x(idx);
            else
                x=x(idx,:);
            end
            D.(vector_fields{j})=x;
        end
    end
end

if isfield(D,'used_normals')
    if isstruct(D.used_normals) && ~isempty(D.used_normals)

        for fl = fieldnames(D.used_normals)'
            D.used_normals.(char(fl)) = D.used_normals.(char(fl))(idx);
        end
 
    elseif iscell(D.used_normals) && ~isempty(D.used_normals)
        D.used_normals  = D.used_normals(idx);
    else
        warning('Removing field: used_normals') %#ok<WNTAG>
        D = rmfield(D,'used_normals');
    end
        
end

%Update remaining matrix fields

matrix_fields={'dat','affy_calls','affy_call','sdat','lsdat','fdat','prectrls','smooth','cbs','cbs_fixed','sm1','sm2','sm2j','sm3','raw','si','flag','hmm','pre_control_sum','sitab','ref','level','scores'};
if isfield(D,'matrix_fields')
    matrix_fields = unique([matrix_fields D.matrix_fields]);
end
for j=1:length(matrix_fields)
    if isfield(D,matrix_fields{j})
        x=D.(matrix_fields{j});
        if ~isempty(x)
            x=x(:,idx);
        end
        D.(matrix_fields{j})=x;
    end
end


%supdat

if isfield(D,'supdat')
    x = D.supdat;
    if ~isempty(x)
        if ndims(x) == 2
            x = x(:,idx);
        end
        if ndims(x) == 3
            x = x(:,idx,:);
        end
        D.supdat = x;
    end
end


if isfield(D,'mergsupdat')
    x = D.mergsupdat;
    if ~isempty(x)
        if ndims(x) == 2
            x = x(:,idx);
        end
        if ndims(x) == 3
            x = x(:,idx,:);
        end
        D.mergsupdat = x;
    end
end

% reorder ziggurat segments in Qs field
if isfield(D,'Qs')
    if islogical(idx)
        idx = find(idx);
    end
    Qfields = intersect({'amp','del','aod','doa'},fieldnames(D.Qs));
    for k = 1:length(Qfields)
        field = Qfields{k};
        % create sample index map (0 means delete)
        map = zeros(1,max(D.Qs.(field)(:,5)));
        map(idx) = 1:length(idx);
        % map samples and delete
        D.Qs.(field)(:,5) = map(D.Qs.(field)(:,5));
        D.Qs.(field)(D.Qs.(field)(:,5)==0,:) = [];
    end
end

% take care of fields that aren't hard-coded
if exist('origsize2','var')
    
    fields = fieldnames(D);
    
    remfields = setdiff(fields,...
        [matrix_fields vector_fields {'used_normals','sample_fields','supdat',...
                        'mergsupdat','dat','gsesc','gsymb','gacc','id','symb',...
                        'affy_call','adat','gsupdat','sdesc','sscale','Qs'}]);

    for fi = 1:length(remfields)
        fl = remfields{fi};
        if size(D.(char(fl)),2) == origsize2
            D.(char(fl)) = D.(char(fl))(:,idx);
        else
            continue
        end
    end
    
end

