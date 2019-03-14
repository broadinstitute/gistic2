function [M,idx1,idx2]=match_string_sets(set1,set2)
%MATCH_STRING_SETS return a sparse array indicating match locations between set1
%and set2
%
%   [M,IDX1,IDX2]=MATCH_STRING_SETS(SET1,SET2) returns a sparse
%   array M, and the ordered pairs (IDX1(k),IDX2(k)) for which M is non-zero.
%   M is a K-by-J matrix, where J is the number of string arrays in SET1,
%   and K is the number of string arrays in SET2. M(j,k)== 1 iff
%   set1(j)==set2(k).  Inputs SET1 and SET2 are cell arrays of strings whose 
%   elements are the strings to match.  (SET1 and SET2 can also be two 
%   dimensional character arrays for which size(set,1) == number of strings 
%   in set).  
%
%   EXAMPLE:
%
%   animals1 = {'cat','fish','hamster','skunk'}
%   animals2 = {'cat','dog','skunk'}
%   [M,mi,mj] = match_string_sets_fast(animals1,animals2)
%   animals1(mi)
%   animals2(mj)
%

% GISTIC software version 2.0
% Copyright (c) 2011, 2016 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


%% Looped implementation
% if iscell(set1)
%   set1=strvcat(set1);
% end
% 
% if iscell(set2)
%   set2=strvcat(set2);
% end
% 
% 
% %ss12 is the set of unique strings
% %strvcat(set1,set2) index by i12 gives ss12
% %ss12(j12) gives strvcat(set1,set2)
% [ss12,i12,j12] = lunique([ strvcat(set1, set2) ],'rows'); 
% n1=size(set1,1);
% n2=size(set2,1);
% %ss12 is the minimal list of words
% 

%M=sparse(n1,n2);
% for i=1:size(ss12,1)
%   idx=find(j12==i); %locations in strvcat containing the ith word
%     
%   idx1=find(idx<=n1);% locations in idx that correspond to set1 
%   idx2=find(idx>n1);% locations in idx that correspond to set2
%     
%   if ~isempty(idx1) & ~isempty(idx2)
%     M(idx(idx1),idx(idx2)-n1)=1; %array with 1 and zeros where word in set1 matches word in set2
%   end
% end

%% Repmat (faster) implementation

if ischar(set1)
    set1 = cellstr(set1);
end

if ischar(set2)
    set2 = cellstr(set2);
end

if size(set1,1)>size(set1,2)
    set1 = set1';
end
if size(set2,1)>size(set2,2)
    set2 = set2';
end


set1cellarray = repmat(set1,length(set2),1);
set2cellarray = repmat(set2',1,length(set1));

matcharray = cellfun(@strcmp,set1cellarray,set2cellarray);



%%

M = sparse(matcharray);

[idx1,idx2,~] = find(M');


