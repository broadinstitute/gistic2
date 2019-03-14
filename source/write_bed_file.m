function write_bed_file(fname,track_name,color_string,regs,C,cyto)


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

f=fopen(fname,'w');
fprintf(f,'track name=%s color=%s\r\n',track_name,color_string);
ampdel={'Any-AP','Any-DP'};
for k=1:length(regs)
  for i=1:length(regs{k})
    [st,chr,bp]=genomic_location(C,{[regs{k}(i).peak_wide_st regs{k}(i).peak_wide_en]},cyto,1);
    fprintf(f,'chr%s\t%d\t%d\t%s\r\n',num2chromosome(chr),bp{1}(1),bp{1}(2),[ampdel{k} num2str(i)]);
  end
end

fclose(f);
