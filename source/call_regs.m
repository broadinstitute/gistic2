function [calls,vals]=call_regs(CL21,regs,ts,no_call_th)


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

if ~exist('no_call_th','var') || isempty(no_call_th)
  no_call_th=ts;
end

calls={};
vals={};

for k=1:size(regs,1)
  for dodel=0:(size(regs,2)-1)
    if ~isempty(regs{k,dodel+1})
      for i=1:length(regs{k,dodel+1})
        %      t1=ts(k);
        if dodel
          t1=-ts{k}(2); % log2(4/3)-1
          nct1=-no_call_th{k}(2);
          ncv=double(full(CL21.dat(regs{k,dodel+1}(i).peak,:)<nct1));
          v=double(full(CL21.dat(regs{k,dodel+1}(i).peak,:)<t1)); % used smooth instead of dat
          v(v==0 & ncv==1)=NaN;
          calls{k,dodel+1}(i,:)=v;
          vals{k,dodel+1}(i,:)=full(CL21.dat(regs{k,dodel+1}(i).peak,:));
        else
          t1=ts{k}(1);
          nct1=no_call_th{k}(1);        
          ncv=double(full(CL21.dat(regs{k,dodel+1}(i).peak,:)>nct1));
          v=double(full(CL21.dat(regs{k,dodel+1}(i).peak,:)>t1)); % used smooth instead of dat
          v(v==0 & ncv==1)=NaN;
          calls{k,dodel+1}(i,:)=v;
          vals{k,dodel+1}(i,:)=full(CL21.dat(regs{k,dodel+1}(i).peak,:));
        end
      end
    else
      calls{k,dodel+1}=[];
      vals{k,dodel+1}=[];     
    end
  end
end
