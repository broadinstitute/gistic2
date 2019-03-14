function D = clean_gistic_input(D,cnv_file,remove_sex,join_segment_size,center)
%CLEAN_GISTIC_INPUT - prepare D-struct of segmented data for GISTIC
%
%   D = clean_gistic_input(D,CNV_FILE,REMOVE_SEX,JOIN_SEGMENT_SIZE,CENTER)
%
% D is a D-struct containing coy number data. Optionally, CNV_FILE names a
% file with excluded regions, REMOVE_SEX is a flag set if sex chromosomes
% should be excluded, join_segment_size is a minimum number of markers
% allowed for a segment, CENTER chooses the method for centering the data:
% 'median', 'mean' or 'none' (default 'median').

% GISTIC software version 2.0
% Copyright (c) 2011, 2016 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

% default join segment size is none
if ~exist('join_segment_size','var') || isempty(join_segment_size)
    join_segment_size = 0;
end

% default centering is median (for compatibility)
if ~exist('center','var')
    center = 'median';
end

% default islog for backward compatibility
if ~isfield(D,'islog')
    D.islog = true;
end

% remove areas of normal copy number variation
if exist('cnv_file','var') && ~isempty(cnv_file)
    D = remove_cnv(D,cnv_file);
end
% check that there are data remaining to process
if size(D.dat,1) == 0
    throw(MException('snp:clean_gistic_input:AllDataRemoved',...
                     'All input data were removed after germline CNV processing.'));
end

% Remove NaN probes
verbose('Removing NaN probes...',20)
nan_idx=find(any(isnan(D.dat),2));
if ~isempty(nan_idx)
    verbose(['Removing ' num2str(length(nan_idx)) ' markers with NaNs'],20); 
    D=reorder_D_rows(D,setdiff(1:size(D.dat,1),nan_idx));
else
    verbose('No markers with NaNs... ',20);
end  
verbose(['Matrix size ' num2str(size(D.dat)) ],20);
% check that there are data left to process
if size(D.dat,1) == 0
  throw(MException('snp:clean_gistic_input:AllDataRemoved',...
                   'All input data were removed after NaN processing.'));
end
  
% remove X,Y chromosome
if remove_sex
    verbose('Removing sex chromosomes...',20);
    RGI = RefGeneInfo;
    keepers = isAutosomal(RGI,D.chrn);
    if isa(D.dat,'SegArray')
        keepers = SegArray(keepers);
    end
    D=reorder_D_rows(D,keepers);
end

% optionally join small segments
if join_segment_size > 0 
    verbose('Merging small segments...',20);
    D = rmfield_if_exists(D,{'cbs','cbs_rl'});
    D.cbs=D.dat;
    D = smooth_cbs(D,join_segment_size);
    D.dat=D.cbs;
    D = rmfield_if_exists(D,{'cbs','cbs_rl'});
end

% optionally center the data
if strcmpi(center,'median')
    % subtract median
    verbose('Median centering data...',20);
    D.medians = median(D.dat,1);
    if D.islog
        % log units: center on medians
        D.dat = D.dat - repmat(D.medians,size(D.dat,1),1);
    else
        % linear units: scale by medians
        D.dat = (2 * D.dat) ./ repmat(D.medians,size(D.dat,1),1) - 2;
    end
elseif strcmpi(center,'mean')
    % subtract median
    verbose('Mean centering data...',20);
    D.means = mean(D.dat,1);
    if D.islog
        % log units: center on means
        D.dat = D.dat - repmat(D.means,size(D.dat,1),1);
    else
        % linear units: scale by means
        D.dat = (2 * D.dat) ./ repmat(D.means,size(D.dat,1),1) - 2;
    end
elseif ~strcmpi(center,'none')
    throw(MException('snp:clean_gistic_input:invalid_parameter',...
                   '''%s'' is not a valid sample-centering option.'));
end
