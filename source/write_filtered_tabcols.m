function nrows = write_filtered_tabcols(file,rows,varargin)
% WRITE_TABCOLS write a tab DLM data file of filtered columns with headers.
%    NROWS = write_filtered_tabcols(FILE,ROWS,{HEAD1,VECT1,FORSPEC1},...,
%    {HEADN,VECTN,FORSPECN}), creates the file FILE and writes N 
%    columns of data described by the N cell array arguments. ROWS 
%    is a numeric or logical index that filters/reorders
%    each column. HEAD1,...,HEADN are column header text strings;
%    VECT1,...,VECTN are column vectors of data, numeric or char; 
%    FORSPEC1,...,FORSPECN are optional format specifiers for
%    converting the column to a string using FPRINTF.  
%
%    the return value is the number of data items written.

% GISTIC software version 2.0
% Copyright (c) 2011,2014 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

nrows = 0;

if ischar(file)
  h = fopen(file,'w');
elseif file > 0
  h = file;
end


if h > 0
    ncols = length(varargin); % number of output columns
    header = '';              % header string
    format = '';              % data format string
    data = cell(1,ncols);     % cell array referencing column data

    % process the ROWS argument
    if ~exist('rows','var') || isempty(rows)
        rows = 1:length(varargin{1}{2});
    end

    % process arguments
    for i = 1:length(varargin)
        % delimit columns
        if i > 1
            header = strcat(header,'\t');
            format = strcat(format,'\t');
        end
        header = strcat(header,varargin{i}(1));
        % build format
        if length(varargin{i}) < 3
            if isnumeric(varargin{i}{2})
                try
                    validateattributes(varargin{i}{2},{'numeric'},{'integer'})
                    forspec = '%d';
                catch
                    forspec = '%g';
                end
            else
                forspec = '%s';
            end
        else
            forspec = varargin{i}{3};
        end
        format = strcat(format,forspec);
        % build argument list
        data{i} = varargin{i}{2}(rows);
    end

    % delimit lines
    header = strcat(header,'\n');
    format = strcat(format,'\n');
    % write header
    fprintf(h,char(header));
    % have to loop over rows because of fprintf limitations :(
    args = cell(1,size(data,2));
    for i = 1:length(data{1});
        for j = 1:length(data);
            if iscell(data{j})
                args{j} = data{j}{i};
            else
                args{j} = data{j}(i);
            end
        end
        % write data
        fprintf(h,char(format),args{:});
        nrows = nrows+1;
    end
    if ischar(file)
        fclose(h);
    end
end
