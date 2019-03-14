function args_table=handle_args(dash_types,args)
% HANDLE_ARGS initialize the structure ARGS_TABLE
%
%   ARGS_TABLE=HANDLE_ARGS(DASH_TYPES,ARGS) populates the structure ARGS_TABLE 
%   with fields and values provided in the cell array of strings, ARGS.  
%   The recognized dash types are given in the cell array of string, DASH_TYPES.
%   The structure fields listed in ARGS are preceeded by a '-'.
%
%   Fields of struct ARGS_TABLE are listed in the odd-indexed elements of ARGS,  
%   and the values corresponding the the fields are listed immediately
%   after the field as even indexed elements.
%
%       Example:  (Typical use for a gp module)
%           a = handle_args({'b','i','o'},varargin}

% GISTIC software version 2.0
% Copyright (c) 2011, 2016 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, Gordon Saksena
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


sargs=strvcat(args);
if isempty(sargs)
    args_table=[];
    for i=1:length(dash_types)
        args_table=setfield(args_table,dash_types{i},[]);
    end
    return
end
dashs=find(sargs(:,1)=='-');            
dashs(end+1)=length(args)+1;            
args_table=[];
for i=1:length(dash_types)
  idx=strmatch(['-' dash_types{i}],args(dashs(1:(end-1))),'exact');
  if ~isempty(idx)
    if dashs(idx+1)==dashs(idx)+1 % no argument
      args_table=setfield(args_table,dash_types{i},[]);
      verbose([dash_types{i} ': []'],20);
    elseif dashs(idx+1)-1==dashs(idx)+1
      args_table=setfield(args_table,dash_types{i},args{dashs(idx)+1});
      verbose([dash_types{i} ': ' regexprep(args{dashs(idx)+1},'\\','\\\\')],20);
    else
      args_table=setfield(args_table,dash_types{i},args((dashs(idx)+1):(dashs(idx+1)-1))); 
      args((dashs(idx)+1):(dashs(idx+1)-1))=regexprep(args((dashs(idx)+1):(dashs(idx+1)-1)),'\\','\\\\');
      verbose([dash_types{i} ': ' sprintf('%s ',args{(dashs(idx)+1):(dashs(idx+1)-1)})],20);      
    end
  else
    args_table=setfield(args_table,dash_types{i},[]);
    verbose([dash_types{i} ': []'],20); 
  end
end
