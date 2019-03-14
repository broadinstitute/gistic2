function C = smooth_cbs(C,n,set_to_val)
%SMOOTH_CBS - join small CN segments to larger ones
%
%  D = SMOOTH_CBS(D,N,SET_TO_VAL)
%
%  The segmented data in D.cbs is processed by joining each segment of
%  length N snps or fewer to the adjacent neighbor closest in copy number. 
%  If SET_TO_VAL is not given or is empty, then the joined segment is a
%  weighted average of the two comprising segments, otherwise, the value of
%  small segments are set to SET_TO_VAL.

% GISTIC software version 2.0
% Copyright (c) 2011, 2016 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

MAXCR = 1e6; % maximum absolute value allowed for copy ratio

if ~isfield(C,'chrn')
    C=add_chrn(C);
end

% default islog
if ~isfield(C,'islog')
    C.islog = true;
end

chrn_rl=runlength(C.chrn);

if ~isfield(C,'cbs_rl') || (isfield(C,'cbs_rl') && isempty(C.cbs_rl))
    % limit copy ratio values to something reasonable
    C.cbs(C.cbs > MAXCR) = MAXCR;
    C.cbs(C.cbs < -MAXCR) = -MAXCR;
    for i=1:size(C.cbs,2)
        C.cbs_rl{i}=runlength(C.cbs(:,i),chrn_rl);
    end
end

verbose('currently not taking care of edge of chromosomes!!!',10);

if size(C.cbs_rl{1},2)<4
    C=add_chrn_to_rl(C);
end

all_chromosomes = min(C.chrn):max(C.chrn);

if ~exist('set_to_val','var') || isempty(set_to_val)
    smooth_chromosomes = all_chromosomes(hist(C.chrn,min(C.chrn):max(C.chrn))>n);
    dontsmooth = setdiff(all_chromosomes,smooth_chromosomes);
    if ~isempty(dontsmooth)
        warning('Number of markers is less than smooth length for the following chromosomes: %s. These chromosomes will not be smoothed.  (Smooth length = %d.)',...
            num2str(dontsmooth),n) %#ok<WNTAG>
    end
else
    smooth_chromosomes = all_chromosomes;
end

%% loop over samples
for i=1:size(C.cbs,2)
    verbose('Smoothing segments for sample %d:',30,i);
    % loop over relevant chromosomes
    for c=smooth_chromosomes
        rl=C.cbs_rl{i};
        in_chr=find(rl(:,4)==c);


        rl_chr=rl(in_chr,:);
        sz=rl_chr(:,2)-rl_chr(:,1)+1;
        small=sz<=n;
        %     large=sz>=large_n;
        if any(small)
            if exist('set_to_val','var') && ~isempty(set_to_val)
                smalls=find(sz<=n & rl_chr(:,3)~=set_to_val);
                while ~isempty(smalls)
                    sj=smalls(1);
                    range=rl_chr(sj,1):rl_chr(sj,2);
                    verbose(['Changing ' genomic_location(C,{range}) ' from ' num2str(rl_chr(sj,3)) ...
                        ' to ' num2str(set_to_val)],50);
                    rl_chr=merge_rl_segs(rl_chr,sj,set_to_val);
                    sz=rl_chr(:,2)-rl_chr(:,1)+1;
                    smalls=find(sz<=n & rl_chr(:,3)~=set_to_val);
                end
                C.cbs_rl{i}=[ rl(1:(in_chr(1)-1),:); rl_chr; rl((in_chr(end)+1):end,:)];
            else
                smalls=find(sz<=n);
                while ~isempty(smalls)
                    rlx=[0 0 Inf 0; rl_chr ; 0 0 Inf 0];
                    delta=[ rlx(smalls+1,3)-rlx(smalls,3) rlx(smalls+2,3)-rlx(smalls+1,3)];
                    [mn,mi]=min(abs(delta),[],1);
                    if mn(1)<mn(2) || isnan(mn(1))  %quick fix to deal with NaNs -- make better
                        C=combine_segs(C,i,in_chr(smalls(mi(1))),in_chr(smalls(mi(1))-1));
                    else
                        C=combine_segs(C,i,in_chr(smalls(mi(2))),in_chr(smalls(mi(2))+1));
                    end
                    rl=C.cbs_rl{i};
                    in_chr=find(rl(:,4)==c);
                    rl_chr=rl(in_chr,:);
                    sz=rl_chr(:,2)-rl_chr(:,1)+1;
                    smalls=find(sz<=n);
                end
            end
        end % if there are small
    end % loop over chromosomes

    % reconstruct array data from segments
    if isa(C.cbs,'SegArray')
        tmp = SegArray.fromRunlength(C.cbs_rl{i});
        C.cbs(:,i) = tmp(SegArray.fromSegments(1,size(C.cbs,1)));
    else
        tmp=derunlength(C.cbs_rl{i});
        C.cbs(:,i) = tmp(1:size(C.cbs,1));
    end
    if size(C.cbs,1)~=size(tmp,1)
        verbose('Assuming removed X chromosome',10);
    end
end % loop over samples


function C=combine_segs(C,i,small_seg,large_seg)
% C is the copy number data structure
% i is the sample index (selects ryn length in C.cbs_rl)
% small_seg is the index of the segment being replaced
% large_seg is the index of the segment it is being joined with
curseg=small_seg;
rl=C.cbs_rl{i};
sz=rl(:,2)-rl(:,1)+1;
curval=rl(curseg,3);
if small_seg>large_seg
    prevseg=large_seg;
    prevval=rl(prevseg,3);
    both=[rl(prevseg,1):rl(prevseg,2),rl(curseg,1):rl(curseg,2)];
    if C.islog
        % if log, use weighted arithmetic mean
        weighted_avg = (prevval*sz(prevseg) + curval*sz(curseg)) ./ (length(both));
    else
        % if CN, use weighted geometric mean
        weighted_avg = prevval^(sz(prevseg)/length(both)) * curval^(sz(curseg)/length(both));
    end
    verbose(['  Joining [' ...
            genomic_location(C,{rl(prevseg,1):rl(prevseg,2)}) ...
            ',' num2str(prevval) '] to <' ...
            genomic_location(C,{rl(curseg,1):rl(curseg,2)}) ...
            ',' num2str(curval) '> --> [' ...
            genomic_location(C,{both}) ...
            ',' num2str(weighted_avg) ']'],50);
    rl=merge_rl_segs(rl,[prevseg curseg],weighted_avg);
    C.cbs_rl{i}=rl;
else
    nextseg=large_seg;
    nextval=rl(nextseg,3);
    both=[rl(curseg,1):rl(curseg,2) rl(nextseg,1):rl(nextseg,2)];
    if C.islog
        % if log, use weighted arithmetic mean
        weighted_avg = (nextval*sz(nextseg) + curval*sz(curseg)) ./ (length(both));
    else
        % if CN, use weighted geometric mean
        weighted_avg = nextval^(sz(nextseg)/length(both)) * curval^(sz(curseg)/length(both));
    end
    verbose(['  Joining <' ...
            genomic_location(C,{rl(curseg,1):rl(curseg,2)}) ...
            ',' num2str(curval) '> to [' ...
            genomic_location(C,{rl(nextseg,1):rl(nextseg,2)}) ...
            ',' num2str(nextval) '] --> [' ...
            genomic_location(C,{both}) ...
            ',' num2str(weighted_avg) ']'],50);
    rl=merge_rl_segs(rl,[curseg nextseg],weighted_avg);
    C.cbs_rl{i}=rl;
end
