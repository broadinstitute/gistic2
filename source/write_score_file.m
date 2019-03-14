function write_score_file(fname,C,p,q,ads,ts)
% WRITE_SCORE_FILE writes an IGV compatible GISTIC scores file
%
%   write_score_file(FNAME,D,P,Q,ADS,TS)
%
% FNAME - output file
% D - copynumber structure (for marker, chr, pos)
% P - UNUSED
% Q - length 2 cell array containing vectors of q-values for amps and dels
%     respectively
% ADS - length 2 cell array containing vectors of G-scores for amps and dels
%       respectively
% TS - length 2 vector of amp and del noise thresholds, default [0.1 0.1]

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


if ~exist('ts','var')
  ts = [0.1 0.1];
end

% open file and write header
f = fopen(fname,'w');
fprintf(f,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',...
    'Type','Chromosome','Start',...
    'End', '-log10(q-value)','G-score', 'average amplitude','frequency');

% write amplification score segments
[amp1,freq1] = ampfreq(C.dat,ts(1));
output_scna_segs(f,'Amp',C,amp1,freq1,q{1},ads{1});

% write deletion score segments
[amp2,freq2] = ampfreq(-C.dat,ts(2));
output_scna_segs(f,'Del',C,amp2,freq2,q{2},ads{2});

% close file
fclose(f);

% SUBFUNCTION: output either amp or del
function output_scna_segs(f,type,C,amp,freq,q,ads)

    % segment the g-score and impose that segmentation on the
    % q-values and chromosomes
    adseg = runlength(ads,C.chrn);
    q1 = q(adseg(:,1));
    chr = num2str(C.chrn(adseg(:,1)));

    % get base positions of marker starts and ends
    st = C.pos(adseg(:,1));
    en = C.pos(adseg(:,2));

    % segmented versions of average frequency and amplitude
    freq_f = full(freq(adseg(:,1)));
    amp_f = full(amp(adseg(:,1)));

    chr_t=char(chr);

    % take the log of the q-value and deal with precision underflow
    logq = -log10(q1);
    logq(isinf(logq)) = nan;
    logq(isnan(logq)) = nanmax(logq);

    % output SCNA segments 
    for i=1:size(adseg,1)
        fprintf(f, '%s\t%s\t%d\t%d\t%f\t%f\t%f\t%f\n', ...
            type, chr_t(i,:), st(i), en(i), logq(i), ...
            adseg(i,3), amp_f(i),freq_f(i));
    end

% SUBFUNCTION: get mean amp level and frequency from the data
% (saves memory in main code)
function [avg_amp,freq] = ampfreq(dat,thresh)
    dat(dat<thresh) = 0;
    freq = mean(dat>0,2);
    freq(freq==0) = eps;
    avg_amp = mean(dat,2) ./ freq;


