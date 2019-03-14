function S = impose_default_value(S,field,value,acceptable_values)
%
% impose_default_value(S,field,value)
%
% S is a struct
% if field does not exist, it is created with the given default value.
% if field does exist but is empty, it is given the default value.
% if field does exist but is nonempty, it is left alone.
%
% Mike Lawrence 2008-06-26
% modified 2009-07-02 to handle required fields
% modified 2011-11-10 to handle accepable values

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


if ~isfield(S,field) || isempty(getfield(S,field))
  if ischar(value) & strcmp(value,'*required*')
    error('%s is a required field of P',field);
  else
    try
      S=setfield(S,field,value);
    catch me
      fprintf('Error setting field "%s"\n',field);
      disp(me);disp(me.message);
    end
  end
end

if exist('acceptable_values','var')
  av = acceptable_values;
  v = getfield(S,field);
  if ischar(v) && isnumeric(av)
    v = str2double(v);
    S = setfield(S,field,v);
  end
  if ~ischar(v) && ischar(av)
    error('%s is assigned value of different type than acceptable_values',field);
  end
  try
    ism = ismember(v,av);
  catch me
    error('%s is assigned value of different type than acceptable_values',field);
  end
  if ~ism
    fprintf('Acceptable values for %s:\n',field); disp(av);
    fprintf('Attempted to set to:\n'); disp(v)
    error('Invalid setting of %s',field);
  end
end
