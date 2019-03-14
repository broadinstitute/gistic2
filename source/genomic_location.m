function [st,chr,bp]=genomic_location(X,i,cyto,enlarge_to_next_snp,use_commas)
% GENOMIC_LOCATION generate text for genomic location(s)
%
%   [ST,CHR,BP]=genomic_location(D,I,CYTO,ENLARGE_TO_NEXT_SNP,USE_COMMAS)
%
% Returns a string in ST representing the genomic location in D indexed by
% I of the form 'chrX:<basepos>(<snp>)' where X is the chromosome,
% <basepos> is the base position, and <snp> is the marker index. The
% chromosome number is returned in CHR and numerical base position in BP.
%
% If I is a cell array containing scalar indices, then a cell array of
% strings is returned in ST, CHR is a row vector of chromosomes, and BP is
% a cell array of scalar locations. If the contents of a cell of I is a
% range of indices, then the corresponding string has the form 
% 'chrX:<startpos>-<endpos>(<startsnp>-<endsnp>)' and the corresponding
% cell in BP will contain the range as a length 2 row vector.
%
% CYTO is an optional argument that must be provided if ENLARGE_TO_NEXT_SNP
% is provided and set. If ENLARGE_TO_NEXT_SNP is set, then the string takes
% the form chrX:<startpos>, 
%! TODOC finish documenting this function

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


%% process optional arguments
if ~exist('use_commas','var') || isempty(use_commas)
  use_commas=0;
end
if ~exist('enlarge_to_next_snp','var') || isempty(enlarge_to_next_snp)
  enlarge_to_next_snp=0;
end

% figure out if Mb units from data
if ~is_mb(X)
    Mb=1;
else
    Mb=1e6;
end

%% process each index
for j=1:length(i)
  if ~iscell(i) && ~enlarge_to_next_snp
    i1=clip_to_range(i(j),[ 1 length(X.chrn)]);
    chr(j)=X.chrn(i1);
    bp(j)=round(X.pos(i1)*Mb);
    st{j}=['chr' num2chromosome(chr(j)) ':' hum_num2str(bp(j),use_commas) '(' num2str(i) ')' ];
  else
    if ~iscell(i)
      snps=[i(j) i(j)];
    else
      snps=i{j};
    end
    i1=clip_to_range(snps,[ 1 length(X.chrn)]);
    chr(j)=X.chrn(min(i1));
    if X.chrn(min(i1))~=X.chrn(max(i1))
      error('crosses chromosomes');
    end

    if enlarge_to_next_snp
      max_in_chr=cyto(max(find(cat(1,cyto.chrn)==X.chrn(min(i1))))).end;
      if max(i1)==length(X.chrn) || X.chrn(max(i1)+1)~=X.chrn(max(i1))
        region_end=max_in_chr;
      else
        region_end=X.pos(max(i1)+1)-1;
      end
      if min(i1)==1 || X.chrn(min(i1)-1)~=X.chrn(min(i1))
        region_start=1;
      else
        region_start=X.pos(min(i1)-1)+1;
      end
      bp{j}=round([region_start region_end]*Mb);
      nbp{j}=round([X.pos(min(i1)) X.pos(max(i1))]*Mb);
    else
      bp{j}=round([X.pos(min(i1)) X.pos(max(i1))]*Mb);
      nbp{j}=bp{j};
    end

    if enlarge_to_next_snp
      st{j}=['chr' num2chromosome(chr(j)) ':' hum_num2str(bp{j}(1),use_commas) '-' hum_num2str(bp{j}(2),use_commas) '(' num2str(min(i1)) ':' ...
             num2str(max(i1)) ',' ...
             'chr' num2chromosome(chr(j)) ':' hum_num2str(nbp{j}(1),use_commas) '-' hum_num2str(nbp{j}(2),use_commas) ...
             ')']; 
    else
      st{j}=['chr' num2chromosome(chr(j)) ':' hum_num2str(bp{j}(1),use_commas) '-' hum_num2str(bp{j}(2),use_commas) '(' num2str(min(i1)) ':' ...
             num2str(max(i1)) ')'];
    end
  end
end

if length(i)==1;
  st=st{1};
end

