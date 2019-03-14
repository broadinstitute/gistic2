function [arm_medians,arm_names,arm_marks] = arm_medians(D,min_arm_marks)
%ARM_MEDIANS - return medians of chromosome arms from a D struct
%
%   [ARM_MEDIANS,ARM_NAMES,ARM_MARKS] = arm_medians(D,MIN_ARM_MARKS)
%
% D is a D struct to which cytoband information has been added.
% MIN_ARM_MARKS specifies the minimum number of markers an arm must have to
% be included. Returns armxsample matrix of medians in ARM_MEDIANS, 
% optionally a cell array of arm names in ARM_NAMES, optionally counts of
% markers

% GISTIC software version 2.0
% Copyright (c) 2011,2014 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

if ~isfield(D,'armn')
    error('D-struct in first argument needs to have cytoband information');
end
if ~exist('min_arm_marks','var') || isempty(min_arm_marks)
    min_arm_marks = 100;
end
Nchr = length(unique(D.chrn));

D.chrn = SegArray(D.chrn);
D.armn = SegArray(D.armn);
arm_medians = nan(Nchr*2,size(D.dat,2));
arm_marks = zeros(Nchr*2,1);
arm_names = repmat({''},Nchr*2,1);
pq = 'pq';
i = 1;

for c=as_row(unique(D.chrn))
    for a=1:2
        picker = D.chrn==c & D.armn==a;
        if any(picker)
            arm_medians(i,:) = nanmedian(D.dat(picker,:));
        end
        arm_marks(i) = sum(picker);
        fprintf('arm %d: %s%s %d markers\n',i,num2chromosome(c),pq(a),sum(picker));
        arm_names{i} = [num2chromosome(c),pq(a)];
        i = i+1;
    end
end
too_small = arm_marks < min_arm_marks;
arm_medians(too_small,:) = [];
arm_names(too_small) = [];
arm_marks(too_small) = [];
