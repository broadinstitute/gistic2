function [Ba Bd] = prepare_B(B)
  
  %% prepares Amp and Del Bs (Ba and Bd, respectively) for
  %atomic_ziggurat_deconstruction

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

    
   if ~isempty(B)
     Ba = B;
     Bd = B;
     
     Ba(:,4) = B(:,4).*(B(:,4)>0);
     Bd(:,4) = -1*(B(:,4).*(B(:,4)<0));
     
     Ba = merge_adj_segs(Ba);
     Bd = merge_adj_segs(Bd);
   else
     Ba = [];
     Bd = [];
   end
