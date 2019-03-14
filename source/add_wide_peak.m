function [s,e,sample_s,sample_e]=add_wide_peak(z,chr_zero,has_unique_rows)

% has_unique_rows also assumes no NaNs

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

if ~exist('has_unique_rows','var')
  has_unique_rows=0; 
end


if ~has_unique_rows
  if nnz(isnan(z))>0
    error('has nan');
  end
  zorig=z;
  zi=[1; 1+find(any(diff(zorig,1,1),2))];
  zj=cumsum([1; any(diff(zorig,1,1),2)]);
  z=z(zi,:);
end

zsc=nanmean(z,2);
[mx,mi]=max(zsc);
rg=find(zsc==mx);
st=min(rg);
en=max(rg);
samples=find(z(mi,:));
zc=z(:,setdiff(1:size(z,2),samples));
zv=z(:,samples);
zsc_c=sum(zc,2);
zsctot=zsc_c+sum(zv,2);

rl=runlength(zv);
if ~iscell(rl)
  rl={rl};
end
rli=find_rl(rl,mi);

top_s=[]; % all samples that have segments that start at the top of peak
bot_s=[]; % all samples that have segments that end at the end of the peak
for i=1:length(rl)
  if rl{i}(rli{i},1)==st
    top_s=[ top_s i];
  end
  if rl{i}(rli{i},2)==en
    bot_s=[ bot_s i];
  end
end


sts=st;
for k=1:length(top_s)
  %zsc=nanmean(z(:,setdiff(1:size(z,2),top_s(k))),2);
  zsc=zsctot-zv(:,top_s(k));
  [mx,mi]=max(zsc);
  rg=find(zsc==mx);
  sts(end+1)=min(rg);
end

ens=en;
for k=1:length(bot_s)
%  zsc=nanmean(z(:,setdiff(1:size(z,2),bot_s(k))),2);
  zsc=zsctot-zv(:,bot_s(k));
  [mx,mi]=max(zsc);
  rg=find(zsc==mx);
  ens(end+1)=max(rg);
end

[spos,sample_s]=min(sts);
if sample_s>1
  sample_s=top_s(sample_s-1);
else
  sample_s=[];
end

[epos,sample_e]=max(ens); %% Changed from min(ens) to max(ens)...source of bug - C.M. 4-25-08
if sample_e>1
  sample_e=bot_s(sample_e-1);
else
  sample_e=[];
end

if ~has_unique_rows
  s=min(find(zj==spos))+chr_zero;
  e=max(find(zj==epos))+chr_zero;
else
  s=spos+chr_zero;
  e=epos+chr_zero;
end
