function Z = atomic_zigg_deconstruction(Bt)
  %% atomic_zigg_deconstruction performs the atomic ziggurat
  %deconstruction operations on a segment dataset Bt, representing the
  %segments for a given chromosome.

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

    
  % Each entry of B represents a copy number segment in a sample.
  % The columns of B are: 1)chrn, 2)st, 3)en, 4)amp, 5)sample, and 6)fract
  % The Fields of Z are: 1)chrn, 2)st, 3)en, 4)amp, 5)sample, 6)cn_st,
  % 7)cn_en, and 8) fract
  % Note that in Z, amp = cn_en - cxn_st
    
    if ~isempty(Bt)
      % Initialize variables
      sample = Bt(1,5);
      chrn = zeros(size(Bt,1)+1,1);
      st = chrn;
      en = chrn;
      fract = chrn;
      cn_st = chrn;
      cn_en = chrn;
      c = 1; % counter variable    
      
      while length(Bt(:,4)) > 1 && max(Bt(:,4)) > 0
        [~,mi] = max(Bt(:,4));
        if mi == 1 %% max is leftmost segment on chromosome
          adj_idx = mi+1;
        elseif mi == size(Bt,1) %% max is rightmost segment on chromosome
          adj_idx = mi-1;
        else %% max is between other segments
          if Bt(mi-1,4) >= Bt(mi+1,4)
            adj_idx = mi-1;
          else
            adj_idx = mi+1;
          end
        end
        diff = Bt(mi,4)-Bt(adj_idx,4);
        if diff>0
          chrn(c) = Bt(mi,1);
          st(c) = Bt(mi,2);
          en(c) = Bt(mi,3);
          fract(c) = Bt(mi,6);
          cn_st(c) = Bt(adj_idx,4);
          cn_en(c) = Bt(mi,4);
          c = c+1;
        end
        
        % Merge adjoining rows (and all further rows with same value)
        kk = min(adj_idx,mi);
        Bt(mi,4) = Bt(adj_idx,4);
        
        % Check if adjacent segs have same value
        while kk < size(Bt,1) && Bt(kk,4) == Bt(kk+1,4)
          Bt(kk,3) = Bt(kk+1,3);
          Bt(kk,6) = Bt(kk,6) + Bt(kk+1,6);
          Bt(kk+1,:) = [];
        end
        
      end
      
      % Check to see if broad level is left
      if Bt(1,4) ~= 0
        chrn(c) = Bt(1,1);
        st(c) = Bt(1,2);
        en(c) = Bt(1,3);
        fract(c) = Bt(1,6);
        cn_st(c) = 0;
        cn_en(c) = Bt(1,4);
        c = c+1;
      end
      
      amp = cn_en-cn_st;
      
      Z = [chrn(1:c-1) st(1:c-1) en(1:c-1) amp(1:c-1) repmat(sample,c-1,1) ...
           cn_st(1:c-1) cn_en(1:c-1) fract(1:c-1)];
      
      %% {'chrn','st','en','amp','sample','cn_st','cn_en','fract'});
    else
      Z = [];
    end
    
            
   
