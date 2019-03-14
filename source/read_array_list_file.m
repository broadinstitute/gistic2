function use_arrays=read_array_list_file(fname,force_text,otherfields)
%READ_ARRAY_LIST_FILE read array list file and return use arrays.  
%
%  USE_ARRAYS = READ_ARRAY_LIST_FILE(FNAME,FORCE_TEXT) returns a structure array
%  USE_ARRAYS with fieldnames specified by the column headings in the array
%  list file, FNAME.  FNAME is a tab-delimited file or an excel file (with
%  .xls extension).  The first line of FNAME contains the tab-delimited
%  field labels ('array' and 'file').  The subsequent lines of FNAME
%  contain the tab-delimited field data.  FORCE_TEXT is an optional logical
%  input that, when true, forces READ_ARRAY_LIST_FILE to treat the info
%  file as tab-delimited.  OTHERFIELDS is a cell array of column names that
%  you wish to include in the read.  If strings are given in OTHERFIELDS,
%  the returned structure will have OTHERFIELDS in addition to the standard
%  fields.


% GISTIC software version 2.0
% Copyright (c) 2011, 2016 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.



%% Read .xls (excel) file

if ~exist('force_text','var')
  force_text = 0;
end

if ~exist('otherfields','var')
    otherfields = [];
end

if strcmp(file_ext(fname),'xls') && ~force_text
  [~,~,raw] = xlsread(fname);
  raw = cellfun(@num2str,raw,'UniformOutput',0);  %make sure everything's type char
  
  nancell = cellfun(@isnan,raw,'UniformOutput',0);
  lencell = cellfun(@length,raw);
  nancell(lencell~=1) = {0};
  
  nan_pos = cellfun(@double,nancell);

  raw(nan_pos) = cellstr(repmat('EMPTY',sum(nan_pos),1));
  f{1} = raw(1,:);
  is_xls = 1;
else
  is_xls = 0;
end
 
%% Read Tab-Delimited data

if ~is_xls
  [f,fid] = read_dlm_file(fname,char(9),1);
end

%% Map fields in file to standard field names



%Info file field names  (later on read columns whos heading begins with inc_')
filefields = regexprep(lower(f{1}),'\(.*\)','');

colnames = {'array','dfile'};
colnames = [colnames otherfields];
includes = strmatch('include in',filefields);
colnames = [colnames f{1}(includes)];

[~,m1,m2] = match_string_sets(colnames,filefields);

if isempty(m2)
    error('No match to column labels in array list file.  May need to include "Array" as column header.')
end

fields = regexprep(colnames,'include in ','inc_');


verbose('...Matching fields to Array List File:',30);
verbose([repmat('      ',length(m1),1),char(f{1}{m2}),repmat('  --->  ',length(m1),1),char(fields{m1})],30);


%% Make Data Structure

if ~is_xls
  form=[repmat('%s',1,length(f{1}))];
  
  %read M columns of sample info file into cells of F_dat (cells are Nx1)
  F_dat=textscan(fid,form,'delimiter','\t','emptyValue',NaN);

  %horizontally concatinate M cells of F_dat (new F_dat is N x M cell)
  F_dat=cat(2,F_dat{:});
  
  empty_pos=find(cellfun('isempty',F_dat));
  
  F_dat(empty_pos)=cellstr(repmat('EMPTY',length(empty_pos),1));


else
  F_dat=raw(2:end,:);
end

%s is structure whos fields are the column names of F_dat
use_arrays=cell2struct(F_dat(:,m2),fields(m1),2);

%% Error catching: Make sure all sample names are unique

[~,ULidx,SLidx] = unique({use_arrays.array});

try ULidx == SLidx;
catch
    N = hist(ULidx,max(ULidx));
    notunique = find(N~=1);
    error(['The array column of the array list file is not unique.'...
        '  Please see rows: ' num2str(notunique)])
end
