function D=add_D_field(D,field_type,fields)
%ADD_D_FIELD - add a re-orderable field to a D structure
%
%   D = add_D_field(D, FIELD_TYPE, FIELDS);
%   
% Add custom FIELDS of type FIELD_TYPE to D. The fields are not
% actually created in D, but they are associated with a reordering
% type. FIELDS is a cell array of name strings. The FIELD_TYPE
% string must be one of 'gene', 'sample', or 'matrix'.
% 

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

if ischar(fields)
  fields={fields};
end

switch field_type
 case 'gene'
  if isfield(D,'gene_fields')
    D.gene_fields=unique([D.gene_fields as_row(fields)]);
  else
    D.gene_fields=[as_row(fields)];
  end
 case 'sample'
  if isfield(D,'sample_fields')
    D.sample_fields=unique([D.sample_fields as_row(fields)]);
  else
    D.sample_fields=[as_row(fields)];
  end
 case 'matrix'
  if isfield(D,'matrix_fields')
    D.matrix_fields=unique([D.matrix_fields as_row(fields)]);
  else
    D.matrix_fields=[as_row(fields)];
  end    
 otherwise 
  error('No such field type');
end

