function [d,res]=dist(X,Y,dist_type,transpose_X,transpose_Y)
% d=dist(X,Y,dist_type)
%    calculates the distance between M rows of X and N rows of Y
%    using the distance metric defined by dist_type.
%    returns a M by N distance matrix.
%    supports NaNs.
%
%    X - data points (rows)
%    Y - data points (rows)
%    dist_type - string: 
%      euclidean, cn_euclid, cosine, other - uses Matlab's pdist
% 
% Gaddy Getz
% Cancer Genomics
% The Broad Institute
% gadgetz@broad.mit.edu
%

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


res={};
if ~exist('dist_type','var')
  dist_type.method='euclidean';
end

if ischar(dist_type)
  tmp.method=dist_type;
  dist_type=tmp;
end

if exist('transpose_X','var') && transpose_X
  X=X';
end

if exist('transpose_Y','var') && transpose_Y
  Y=Y';
end

if nnz(isnan(X))+nnz(isnan(Y))>0
  has_nan=1;
else
  has_nan=0;
end

if isempty(Y)
  switch dist_type.method
   case 'euclidean'
    if ~has_nan
      d=dna_dist(X);
    else
      d2=nan_mat_prod(X.^2,~isnan(X'))+nan_mat_prod(~isnan(X),(X.^2)')-2*nan_mat_prod(X,X');
      d2(d2<0)=0;
      ngood=double(~isnan(X))*double(~isnan(X'));
      d2=size(X,2)*d2./ngood;
      d=sqrt(d2);
    end
   case 'nan_corr'
    d=1-nan_corr(X,X);
   case 'cn_euclid'
    if ~has_nan
      d=dist(dna_norm(X,1),[],'euclidean');
    else
      d2=2*(1-nan_corr(X,X));
      d2(d2<0)=0;
      d=sqrt(d2);      
    end
   case 'nanjaccard'
    d=1-jaccard(X,X);
   otherwise
    if isfield(dist_type,'params')
      if iscell(dist_type.params)
        d=squareform(pdist(X,dist_type.method,dist_type.params{:}));
      else
        d=squareform(pdist(X,dist_type.method,dist_type.params));
      end
    else
      d=squareform(pdist(X,dist_type.method));
    end
  end
else
  switch dist_type.method
   case 'euclidean'
    sx2=sum(X.^2,2);
    sy2=sum(Y.^2,2);
    d=repmat(sx2,1,size(Y,1))+repmat(sy2',size(X,1),1)-2*X*Y';
    d(d<0)=0;
    d=sqrt(d);
   case 'euclidean_sqrd'
    sx2=sum(X.^2,2);
    sy2=sum(Y.^2,2);
    d=repmat(sx2,1,size(Y,1))+repmat(sy2',size(X,1),1)-2*X*Y';
    d(d<0)=0;
   case 'euclidean_after_med_subtract'
    X=X-repmat(nanmedian(X,2),1,size(X,2));
    Y=Y-repmat(nanmedian(Y,2),1,size(Y,2));
    sx2=nanmean(X.^2,2)*size(X,2);
    sy2=nanmean(Y.^2,2)*size(Y,2);
    if has_nan
      d=repmat(sx2,1,size(Y,1))+repmat(sy2',size(X,1),1)-2*nan_mat_prod(X,Y')./(double(~isnan(X))*double(~isnan(Y')))*size(X,2);
      d(d<0)=0;
    else
      d=repmat(sx2,1,size(Y,1))+repmat(sy2',size(X,1),1)-2*X*Y';
      d(d<0)=0;
    end
   case 'seuclidean'
    dist_type.method='seuclidean_sqrd';
    d=sqrt(dist(X,Y,dist_type));
   case 'seuclidean_sqrd'
    w=dist_type.inv_sig2;
    sx2=sum(repmat(w,size(X,1),1).*X.^2,2);
    sy2=sum(repmat(w,size(Y,1),1).*Y.^2,2);
    d=repmat(sx2,1,size(Y,1))+repmat(sy2',size(X,1),1)-2*X*diag(w)*Y';
    d(d<0)=0;
   case 'mahalanobis'
    dist_type.method='mahalanobis_sqrd';
    d=sqrt(dist(X,Y,dist_type));
   case 'mahalanobis_sqrd'
    W=dist_type.inv_sig;
    sx2=X'*W*X;
    sy2=Y'*W*Y;
    d=repmat(sx2,1,size(Y,1))+repmat(sy2',size(X,1),1)-2*X*W*Y';
    d(d<0)=0;
   case 'nan_corr'
    d=1-nan_corr(X,Y);
   case 'cn_euclid'
    if ~has_nan
      d=dist(dna_norm(X,1),dna_norm(Y,1),'euclidean');
    else
      d2=2*(1-nan_corr(X,Y));
      d2(d2<0)=0;
      d=sqrt(d2);      
    end
   case 'nanjaccard'
    d=1-jaccard(X,Y);    
   case 'cosine'
    X=X./repmat(sqrt(sum(X.^2,2)),1,size(X,2));
    Y=Y./repmat(sqrt(sum(Y.^2,2)),1,size(Y,2));
    d=1-X*Y';
   case 'fisher'
    X1=double(X==1);
    Y1=double(Y==1);
    if nnz(isnan(X))==0 && nnz(isnan(Y))==0
      tot=size(X1,2)*ones(size(X1,1),size(Y1,1));
      sx=repmat(sum(X1,2),1,size(Y1,1));
      sy=repmat(sum(Y1,2),1,size(X1,1))';
    else
      XnotNaN=double(~isnan(X));
      YnotNaN=double(~isnan(Y));
      tot=XnotNaN*YnotNaN';
      sx=X1*YnotNaN';
      sy=XnotNaN*Y1';
    end
    both=X1*Y1';
    ac=both; bc=sx-both; cc=sy-both; dc=tot-sx-sy+both;
    if isfield(dist_type,'laplace')
      ac=ac+dist_type.laplace;
      bc=bc+dist_type.laplace;
      cc=cc+dist_type.laplace;
      dc=dc+dist_type.laplace;
    end
    [p1,p2]=fisher_exact_test(ac(:),bc(:),cc(:),dc(:));
    res={ac,bc,cc,dc,reshape(p2,size(ac,1),size(ac,2))};
    if isfield(dist_type,'fdr')
      fdr=calc_fdr_value(p2');
      p2(fdr>dist_type.fdr)=1;
    end
    d=-log2(reshape(p2,size(ac,1),size(ac,2)));
    s=double( ac./(ac+bc+eps)>cc./(cc+dc+eps) );
    d=d.*(2*s-1);
   case 'rank'
    m=size(X,1);
    n=size(Y,1);
    d=size(X,2);
    X=X';
    Y=permute(Y,[ 2 3 1]);

    tmp=tril(nan(d,d),0);
    ds=repmat(1:d,d,1);
    pairs=find(~isnan(tmp));
    da=ds(pairs);
    ds=ds';
    db=ds(pairs);

    SX=sign(X);
    SY=sign(Y);

%    FA=repmat(SX(da,:),[1 1 n]).*repmat(SY(da,:,:),[1 m 1]);
%    FB=repmat(SX(db,:),[1 1 n]).*repmat(SY(db,:,:),[1 m 1]);

    SDX=sign(X(da,:)-X(db,:));
    SDY=sign(Y(da,:,:)-Y(db,:,:));
    
%    FD=repmat(SDX,[1 1 n]).*repmat(SDY,[1 m 1]);
%    d=squeeze(sum( (FA==1).*(FB==1)+(FD==1).*((FA==1)+(FB==1))+(FA==0).*(FB==0),1));
    
    SDX=repmat(SDX,[1 1 n]);
    SDY=repmat(SDY,[1 m 1]);
    SXA=repmat(SX(da,:),[1 1 n]);
    SXB=repmat(SX(db,:),[1 1 n]);
    SYA=repmat(SY(da,:,:),[1 m 1]);
    SYB=repmat(SY(db,:,:),[1 m 1]);
    
    good_pairs={[1 1 0; 1 1 -1],[1 1 0; 1 1 0],[1 1 0; 1 1 1],[1 0 1;1 1 1],[1 0 1;1 -1 1],...
                [1 -1 1;1 -1 1],[0 -1 1; 1 -1 1],[0 -1 1;-1 -1 1],...
                [0 0 0; 1 1 -1],[0 0 0; 1 1 0],[0 0 0; 1 1 1],[0 0 0; 0 0 0],[0 0 0; 1 -1 1],...
                [0 0 0; -1 -1 -1],[0 0 0; -1 -1 0],[0 0 0; -1 -1 1],...
                [0 1 -1; 1 1 -1],[0 1 -1;-1 1 -1],[-1 1 -1;-1 1 -1],[-1 0 -1;-1 1 -1],[-1 0 -1;-1 -1 -1],...
                [-1 -1 0;-1 -1 -1],[-1 -1 0; -1 -1 0],[-1 -1 0;-1 -1 1]};

    dst=zeros(m,n);
    for i=1:length(good_pairs)
      t=squeeze(sum( (SXA==good_pairs{i}(1,1)) & (SXB==good_pairs{i}(1,2)) & (SDX==good_pairs{i}(1,3)) & ...
                       (SYA==good_pairs{i}(2,1)) & (SYB==good_pairs{i}(2,2)) & (SDY==good_pairs{i}(2,3)) ,1));
%      good_pairs{i}
%      t
%      pause
      if n==1
        t=t';
      end
      dst=dst+t;
    end
    mx=squeeze(sum(abs(SX),1));
    mx=mx.*(mx-1)/2+(d-mx).*(d-mx-1)/2;
    dst=dst./repmat(mx',1,size(dst,2));
    dst=1-dst;
    
    d=dst;
   otherwise
    d=dist([X; Y],[],dist_type);
    nx=size(X,1);
    d=d(1:nx,nx+(1:size(Y,1)));
  end
end




