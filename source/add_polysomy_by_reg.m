function [D,supids]=add_polysomy_by_reg(D,regs,t,p,no_call_th)


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

if ~exist('p','var')
  p=0.5;
end

if ~exist('t','var')
  t=[0.3 0.3];
end

if ~exist('no_call_th','var')
  no_call_th=[0.3 0.3];
end

if ~isfield(D,'chrn')
  D=add_chrn(D);
end

supids=[];
ampdel={'A','D'};
for k=1:length(regs)
  for i=1:length(regs{k})
    if regs{k}(i).broad
      name=[num2str(i) ampdel{k} 'B'];
      if k==1
        v_no_call=double(mean(D.dat(regs{k}(i).st:regs{k}(i).en,:)>no_call_th(1),1)>p);
        v=double(mean(D.dat(regs{k}(i).st:regs{k}(i).en,:)>t(1),1)>p);
        v(v==0 & v_no_call==1)=NaN;
      else
        v_no_call=double(mean(D.dat(regs{k}(i).st:regs{k}(i).en,:)<-no_call_th(2),1)>p);
        v=double(mean(D.dat(regs{k}(i).st:regs{k}(i).en,:)<-t(2),1)>p);
        v(v==0 & v_no_call==1)=NaN;
      end        
      [D,si]=add_D_sup(D,name,name,v,'cols');
      supids=[supids; si];
    end
  end
end

