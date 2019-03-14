function [log_hd xamp ylen] = generate_2d_hists(QA,QD,xamp,ylen,pseudocount,do_plot)
  

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

  if ~exist('xamp','var') || isempty(xamp)
    xamp = -2:.08:2;
  end
  
  if ~exist('ylen','var') || isempty(ylen)
    ylen = 0:0.04:2;
  end
  
  if ~exist('pseudocount','var') || isempty(ylen)
    pseudocount = 0;
  end
  
  if ~exist('do_plot','var') || isempty(do_plot)
    do_plot = 0;
  end
  
  Q = cat(1,QA,QD);
  % Generate histograms
  hd = hist2d(Q(:,4),Q(:,8),xamp,ylen);
  
  % Add pseudocounts
  hd1 = hd+pseudocount/100*sum(sum(hd));
  
  % Normalize and take log
  hd2 = hd1/sum(sum(hd1));
  log_hd = log(hd2);
  
  if do_plot
    figure()
    bar3(log(hd2)+15);
    set(gca,'XTick',6:5:51);
    set(gca,'YTick',1:10:51);
    set(gca,'XTickLabel',{.2,.4,.6,.8,1,1.2,1.4,1.6,1.8,2})
    set(gca,'YTickLabel',{-2,-1.2,-.4,.4,1.2,2})
    xlabel('Length (fract of chr arm)')
    ylabel('Amplitude')
    zlabel('log(freq)');
    title('Distribution of segments as function of length and amplitude')
  end
  
