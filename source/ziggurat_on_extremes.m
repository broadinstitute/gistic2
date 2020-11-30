function [QAS QDS] = ziggurat_on_extremes(B,min_level,max_level,hd,xamp,ylen)

% GISTIC software version 2.0
% Copyright (c) 2011-2017 Gad Getz, Rameen Beroukhim, Craig Mermel,
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, Gordon Saksena
% All Rights Reserved.
% (See the accompanying LICENSE file for licensing details.)
  
%% Performs ziggurat on the extremes and returns two matrices, QAS
%('amplified' segments) and QDS ('deleted' segments).  These two
%matrices are in turn used to find the 'broad level' for each
%chromosome (or arm) using the iterative_ziggurat function.

  BA = B;
  BD = B;
  
  BA(:,4) = B(:,4)-min_level; %% BA is entirely above zero
  BD(:,4) = -1*(B(:,4)-max_level); %% BD entirely below zero, then inverted
  
  BA = merge_adj_segs(BA);
  BD = merge_adj_segs(BD);
  
  QAS = atomic_zigg_deconstruction(BA);
  QDS = atomic_zigg_deconstruction(BD);
  
  if ~isempty(QAS)
    QAS(:,6:7) = QAS(:,6:7)+repmat(min_level,size(QAS,1),2);
    QAS(:,9) = score_ziggs_by_table(QAS,hd,xamp,ylen);
  end
  
  if ~isempty(QDS)
    QDS(:,6:7) = -1*(QDS(:,6:7)-repmat(max_level,size(QDS,1),2));
    QDS(:,4) = -1*QDS(:,4);
    QDS(:,9) = score_ziggs_by_table(QDS,hd,xamp,ylen);
  end
  
