function [B,E,C,V] = getSegments(S,F)
% return segments from SegArray (inverse of SegArray.fromSegments)
%
%   [B,E,C,V] = getSegments(S,F)
%
% Return Begin, End, Column, and Value for SegArray S. The optional
% filler value F assumes F is the filler and does not return segmnents 
%with value F.

V = getvals(S);
[B,C] = find(getbpts(S));
% construct end from beginning
E = B-1;
E(E==0) = size(S,1);
E = [E(2:end);size(S,1)];

% if "filler" value provided, eliminate segments with that value
if exist('F','var') & ~isempty(F)
    assert(isscalar(F)); %!!! TODO real exception
    keepers = (V ~= F);
    B = B(keepers);
    E = E(keepers);
    C = C(keepers);
    V = V(keepers);
end
