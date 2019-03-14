function [rg]=add_annotation(rg,annot_file,mark)
% ADD_ANNOTATION add highlighted genes to reference genome

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


     % rg.symbol is annotated version of rg.symb 
     for i=1:length(rg)
         rg(i).symbol=rg(i).symb;
     end

     % annotate if a file argument is provided and the file exists
     if ~isempty(annot_file) && exist(annot_file,'file')
         annot_file=read_dlm_file('annot_file');
         % default mark
         if ~exist('mark','var') || isempty(mark)
             mark='**';
         end
         % extract locus -> symbol mapping
         for i=2:length(annot_file)
             fsymb{i-1}=annot_file{i}{1};
             flocusid(i-1)=str2num(annot_file{i}{3});
         end
         % apply mapping
         [Nt,n1,n2]=match_string_sets_hash(cellstr(num2str(flocusid')),cellstr(num2str([rg.locus_id]')));
         missing1a=setdiff(1:length(fsymb),unique(n1));
         un2=unique(n2);
         for i=1:length(un2)
             rg(un2(i)).symbol=[ rg(un2(i)).symbol mark];
         end
     end
end
