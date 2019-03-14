function [Q QA QD] = iterative_ziggurat(QAS,QDS,cur_level,hd,xamp,ylen)


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

if ~isempty(QAS)
    QA = QAS;
    % amplifications must be at minimum the current level
    QA(:,6) = max(QAS(:,6),cur_level);
    changedA = QAS(:,6)-QA(:,6) ~= 0;
    % calculate event amplitude and keep amplifications
    QA(:,4) = QA(:,7)-QA(:,6);
    keep_A = QA(:,4) > 0;
    % rescore the amplifications whose amplitude changed
    changedA = changedA & keep_A;
    if any(changedA)
        QA(changedA,9) = score_ziggs_by_table(QA(changedA,:),hd,xamp,ylen);
    end
    QA = QA(keep_A,:);
else
    QA = [];
end
  
if ~isempty(QDS)
    QD = QDS;
    % deletions must be at maximum the current level
    QD(:,6) = min(QDS(:,6),cur_level);
    changedD = QDS(:,6)-QD(:,6) ~= 0;
    % calculate event amplitude and keep deletions
    QD(:,4) = QD(:,7)-QD(:,6);
    keep_D = QD(:,4) < 0;
    % rescore the deletions whose amplitude changed
    changedD = changedD & keep_D;
    if any(changedD)
        QD(changedD,9) = score_ziggs_by_table(QD(changedD,:),hd,xamp,ylen);
    end
    QD = QD(keep_D,:);
else
    QD = [];
end
% first ouutput is a concatenation of amplification and deletion events
Q = cat(1,QA,QD);
  
  
  
