function [handle y_hat b stats r p] = tufte_plots(x_values,y_values,marker_labels,handle,qA,qD,do_regression)
      

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

  if ~exist('handle','var') || isempty(handle)
    handle = figure();
  end
  
  if ~exist('do_regression','var') || isempty(do_regression)
    do_regression = 1;
  end
  
  if size(x_values,1) < size(x_values,2)
    x_values = x_values';
  end
  
  if size(y_values,1) < size(y_values,2)
    y_values = y_values';
  end
  
  if size(x_values,2) ~= size(y_values,2)
    error('Must have same number of datasets in x and y')
  end
  
  figure(handle);
        
  scatter(x_values,y_values,'Marker','none')
  [r,p] = corr(x_values,y_values);
    
  [b stats] = robustfit(x_values,y_values);
  y_hat = b(2)*x_values+b(1);
  
  if do_regression
    legend(['r = ' num2str(r,3) ' (p = ' num2str(p,3) ')']);
    line(x_values,y_hat,'Color','red');
  end
  
  for i=1:length(marker_labels)
    text(x_values(i),y_values(i),marker_labels{i});
  end
  
  
