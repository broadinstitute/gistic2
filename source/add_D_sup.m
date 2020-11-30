function [D,supidx]=add_D_sup(D,supacc,supdesc,supdat,rc,allowoverwrite)
%ADD_D_SUP Add a line of supplemental data to data structure.
%
%   [D,SUPIDX] = ADD_D_SUP(D,SUPACC,SUPDESC,SUPDAT,RC,ALLOWOVERWRITE)
%   returns a new data structure D with one row of updated supplemental
%   information. SUPACC is the short field (it is appended to .SUPACC or
%   .GSUPACC). SUPDESC is the descriptive field (appended to .SUPDESC or
%   .GSUPDESC). SUPDAT is the supplemental data given as a vector with
%   length equal to the number of samples or genes (appended to .SUPACC or
%   .GSUPACC). RC is the type of supplemental information, either 'col' for
%   sample information or 'row' for gene information.  Setting 'col'
%   appends the information to .SUP* and 'row' appends the data to .GSUP*.
%   ALLOWOVERWRITE has a default value of 1.  If ALLOWOVERWRITE is 1, then
%   if a SUPACC already exists in the .supacc field, the corresponding
%   .supdat will be overwritten with SUPDAT.
%
%   Example: If D is a data structure containing information on 6 samples,
%            and we define the following:
%
%                   acc = 'Gender: 1-M/2-F' desc = 'Gender 1-Male/2-Female'
%                   
%                   data = [1 1 2 2 2 1]
%
%            then [D,supidx] = add_D_sup(D,acc,desc,data,'cols') :
%                  - Appends a row to the char array D.supacc as given by acc
%                  - Appends a row to the char array D.supdesc as given by desc
%                  - Appends a row to the double array D.supdat as given by
%                  data

% GISTIC software version 2.0
% Copyright (c) 2011-2017 Gad Getz, Rameen Beroukhim, Craig Mermel,
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, Gordon Saksena
% All Rights Reserved.
% (See the accompanying LICENSE file for licensing details.)

if ~exist('rc','var')
    rc='cols';
end

if ~exist('allowoverwrite','var')
    allowoverwrite =1;
end

if ~iscell(supacc)
  supacc=cellstr(supacc);
end
supacc=as_column(supacc);

if ~iscell(supdesc)
  supdesc=cellstr(supdesc);
end
supdesc=as_column(supdesc);

if is_col(rc)
  %add sample supplemental info (columns)
  if ~isfield(D,'supacc'), D.supacc=[]; end;  %supacc should be char arrays, not cells!!
  if ~isfield(D,'supdesc'), D.supdesc=[]; end; %see above
  if ~isfield(D,'supdat'), D.supdat=[]; end;    
  
  %quick fix to prevent a duplicate SUPACC from being added
  
  if allowoverwrite
    supidx = find_supid(D,supacc);
    if ~isempty(supidx)
      if ndims(D.supdat) == 2
        D.supdat(supidx,:) = supdat;
      end
      if ndims(D.supdat) == 3
        if ndims(supdat) == 2
          supdat = repmat(supdat,[1,1,size(D.supdat,3)]);
        end
        
        D.supdat(supidx,:,:) = supdat;
      end
      return
    end
    
  end

  %%% Add supacc
  if iscell(D.supacc)
    D.supacc=[ D.supacc; supacc];
  else
    D.supacc=strvcat(D.supacc,strvcat(supacc));
  end
  
  %% Add supdesc
  
  if isempty(supdesc) && ~isempty(supacc)
    supdesc=supacc;
  end
  if iscell(D.supdesc)
    D.supdesc=[ D.supdesc; supdesc];
  else
    D.supdesc=strvcat(D.supdesc,strvcat(supdesc));
  end
  
  
  %% Add supdat
  
  
  if ndims(D.supdat) == 3
    if ndims(supdat) == 2
      supdat = repmat(supdat,[1,1,size(D.supdat,3)]);
    end
  end
  
  D.supdat=[D.supdat; supdat];
  if isempty(supidx)
    supidx=(size(D.supdat,1)-size(supdat,1)+1):size(D.supdat,1);
  end
  
else
  %add gene supplemental info (rows)
  if ~isfield(D,'gsupacc'), D.gsupacc={}; end;
  if ~isfield(D,'gsupdesc'), D.gsupdesc={}; end;
  if ~isfield(D,'gsupdat'), D.gsupdat=[]; end;
  
  if allowoverwrite
    supidx = find_supid(D,supacc,'rows');
    if ~isempty(supidx)
      if ndims(D.gsupdat) == 2
        D.gsupdat(supidx,:) = supdat;
      end
      if ndims(D.gsupdat) == 3
        if ndims(supdat) == 2
          supdat = repmat(supdat,[1,1,size(D.gsupdat,3)]);
        end
        
        D.gsupdat(supidx,:,:) = supdat;
      end
      return
    end
    
  end
  
  if iscell(D.gsupacc)
    D.gsupacc=[ D.gsupacc; supacc];
  else
    D.gsupacc=strvcat(D.gsupacc,strvcat(supacc));
  end
  
  if iscell(D.gsupdesc)
    D.gsupdesc=[ D.gsupdesc; supdesc];
  else
    D.gsupdesc=strvcat(D.gsupdesc,strvcat(supdesc));
  end
  
  if ndims(D.gsupdat) == 3
    if ndims(supdat) == 2
      supdat = repmat(supdat,[1,1,size(D.gsupdat,3)]);
    end
  end
  
  D.gsupdat=[D.gsupdat; supdat];
  
  if isempty(supidx)
    supidx=(size(D.gsupdat,1)-size(supdat,1)+1):size(D.gsupdat,1);
  end
end

