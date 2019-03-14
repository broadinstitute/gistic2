function [h,gr]=subplotgrid(gr,i,j,change_subplot)


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

mini=min(i);
minj=min(j);
if isnan(gr{mini,minj}.handle)
  maxi=max(i);
  maxj=max(j);
  posmin=gr{mini,minj}.position;
  posmax=gr{maxi,maxj}.position;
  pos=[ posmin(1) posmax(2) posmax(1)+posmax(3)-posmin(1) posmin(2)+posmin(4)-posmax(2)];
  gr{mini,minj}.handle=subplot('position',pos);
  h=gr{mini,minj}.handle;
  axis off
  box off
else
  if exist('change_subplot','var') && ~change_subplot
    c=gca;
    h=subplot(gr{mini,minj}.handle);
    subplot(c);
  else
    h=subplot(gr{mini,minj}.handle);
  end    
end
