function gistic_broad_focal_heatmaps(Dpath,params,broad_out_fname,focal_out_fname,filetype)
%GISTIC_BROAD_FOCAL_HEATMAPS create heatmaps for reconstructed broad/focal genomes
%
%   gistic_broad_focal_heatmaps(DPATH,PARAMS,BROAD_OUT_FNAME,FOCAL_OUT_FNAME,FILETYPE);
%
% Create heatmaps of broad and focal genomes from a GISTIC D-structure and
% save them as files. DPATH is the file path to the copy number D-struct 
% saved by GISTIC. BROAD_OUT_FNAME and FOCAL_OUT_FNAME are paths to the
% files to save the broad and focal heatmap files to. FILETYPE is the
% matlab image file type identifier, usually the same as the extension. The
% default FILETYPE is 'png'
%
% PARAMS is the same parameter structure passed to GISTIC, the specific fields used are:
%
%   PARAMS.t_amp and PARAMS.t_del - amplification and deletion noise
%     thresholds
%   PARAMS.broad_len_cutoff - broad length event cutoff threshold, in
%     chromosome arm units
%

% default file type in png
if ~exist('filetype','var') || isempty(filetype)
    filetype = 'png';
end

% load the saved GISTIC data, including the events in the Qs field
D = load_D(Dpath);
[nmarkers nsamples] = size(D.dat);

%% broad heatmap

% reconstruct focal genome
recon_params = struct('broad_or_focal','focal',...
                      'broad_len_cutoff',params.broad_len_cutoff,...
                      't_amp',params.t_amp,...
                      't_del',params.t_del,...
                      'column_to_add',4,...
                      'use_segarray',true,...
                      'rows',nmarkers,...
                      'cols',nsamples);

focals = reconstruct_genomes(D.Qs,recon_params);
% replace data in local D (too save memory) with focal data
D.dat = focals.amp - focals.del;
figure;
display_D_heatmap(D);
% save focal heatmap image
saveas(gcf,focal_out_fname,filetype);

%% focal heatmap

% reconstruct broad_genome
recon_params.broad_or_focal = 'broad';
broads = reconstruct_genomes(D.Qs,recon_params);

% replace data in local D (too save memory) with broad data
D.dat = broads.amp - broads.del;
figure;
display_D_heatmap(D);
% save broad heatmap image
saveas(gcf,broad_out_fname,filetype);


