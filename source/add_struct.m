function st=add_struct(a,b)


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

if isempty(a)
  st=b;
elseif prod(size(a))==1
  st=add_struct_single(a,b);
elseif any(size(a)-size(b))
  error('a and b are not in the same size');
else
  for i=1:size(a,1)
    for j=1:size(a,2)
      st(i,j)=add_struct_single(a(i,j),b(i,j));
    end
  end
end

function st=add_struct_single(a,b)

st=a;
fn=fieldnames(b);
for i=1:size(fn,1)
  if isfield(st,fn{i})
    if isstruct(getfield(st,fn{i}))
      %	disp(['merge ' fn{i}])
      st=setfield(st,fn{i},add_struct(getfield(st,fn{i}),...
                                      getfield(b,fn{i})));
    elseif iscell(getfield(st,fn{i}))
      %	disp(['merge ' fn{i}])
%      warning('new feature w/cell....check if works correctly');
      st=setfield(st,fn{i},[ getfield(st,fn{i}) getfield(b,fn{i})]);
    else
      %	disp(['overwrite ' fn{i}])
      st=setfield(st,fn{i}, getfield(b,fn{i}));
    end
  else
    %      disp(['add ' fn{i}])
    st=setfield(st,fn{i},getfield(b,fn{i}));
  end
end
