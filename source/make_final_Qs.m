function [QA QD QAOD QDOA] = make_final_Qs(Q)
  

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

  QA = Q(find(Q(:,4) > 0 & Q(:,6) >= 0),:);
  QD = Q(find(Q(:,4) < 0 & Q(:,6) <= 0),:);
  
  QAOD = Q(find(Q(:,4) > 0 & Q(:,7)<=0),:);
  QDOA = Q(find(Q(:,4) < 0 & Q(:,7)>=0),:);
  
  spans_zeroA = find(Q(:,4) >0 & Q(:,6).*Q(:,7) < 0);
  spans_zeroD = find(Q(:,4) < 0 & Q(:,6).*Q(:,7)<0);
  
  new_a = cell(1,length(spans_zeroA));
  new_aod = cell(1,length(spans_zeroA));
  for i=1:length(spans_zeroA)
    cur_seg = Q(spans_zeroA(i),:);
    aod_seg = [cur_seg(1:6) 0 cur_seg(8:10)];
    amp_seg = [cur_seg(1:5) 0 cur_seg(7:10)];
    aod_seg(4) = aod_seg(7)-aod_seg(6);
    amp_seg(4) = amp_seg(7)-amp_seg(6);
    new_a{i} = amp_seg;
    new_aod{i} = aod_seg;
  end
  
  new_d = cell(1,length(spans_zeroD));
  new_doa = cell(1,length(spans_zeroD));
  for i=1:length(spans_zeroD)
    cur_seg = Q(spans_zeroD(i),:);
    doa_seg = [cur_seg(1:6) 0 cur_seg(8:10)];
    del_seg = [cur_seg(1:5) 0 cur_seg(7:10)];
    doa_seg(4) = doa_seg(7)-doa_seg(6);
    del_seg(4) = del_seg(7)-del_seg(6);
    new_d{i} = del_seg;
    new_doa{i} = doa_seg;
  end
  
  QA = cat(1,QA,new_a{:});
  QD = cat(1,QD,new_d{:});
  
  QAOD = cat(1,QAOD,new_aod{:});
  QDOA = cat(1,QDOA,new_doa{:});
  
  
