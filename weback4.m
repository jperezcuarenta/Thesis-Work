% This code pertains to Figures 3.1, 3.2, 3.3, and 3.4 in thesis manuscript

% Read warning before running script

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: 
% f: text for which image to consider 
% e.g. 'Squares', 'Circles', or 'UserSelect'
% if third option then open only .png images
% alpha: which scales to favor

% Output:
% Four figures,
% Fig1: Image being studied 
% Fig2: Average of Sf curves
% Fig3: Region Wise Classification
% Fig4: K-means classification
% Fig5: Agglo classification 

% Notes:
% 1) Can use hierarchichal if there is no intuition 
% on which scales are of interest 
% 2) Want to test diff between k means and 
% taking pixels corresponding to regions in
% average Sf
% 3) Required scripts:
%       yellowSquare.m, yellowCircles.m, weback5.m
% 4) Parameters used in thesis figures: 
%       Region: [65 194]
%       txt = 'Squares';
%       alphaIn = 0;
% 5) Required .mat files:
%       'kMean_Copy.mat','kMeansMatrix.mat','kMeansMatrix_3Clusters.mat'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WARNING:
% DO NOT USE NEW GENERATED MATRICES
% USE SAVED ONES IN WORKSPACE 
% tScaleMatrix, kMeansMatrix, AggoMatrix

% Use kMean_Copy (load from .mat) 
% to Produce k means clustering graph for thesis 
% since I have already changed the colors to match the other figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(txtForImage,'Squares')
% square example 
f = yellowSquare(1);
elseif strcmp(txtForImage,'Circles')
% circle example 
f = yellowCircles(1);
elseif strcmp(txtForImage,'UserSelect') 
% images from folder 
 folderName = 'C:\Users\Jesus Perez Cuarenta\Documents\MATLAB\Thesis\Input Images';
 ext='*.png' ; 
 [filename, path] = uigetfile(fullfile(folderName,ext)) ;
 ImageData = imread(fullfile(path, filename));
 f = imresize(ImageData, [256 256]);
% f = rgb2gray(f);
 f = im2double(f);
end

% pre process data 
f = imadjust(f,stretchlim(f),[]);

figure; imshow(f); 
[M, N] = size(f);

%%%%%%%%%%%%%%%%%%%%%%%%
% not needed if choosing alpha manually 
% alphaIn = weback5(f)
%%%%%%%%%%%%%%%%%%%%%%%%

% Compute scales
Sf = localScale(f,alphaIn);
n = 255;
tauAxis = 0:n;
% Compute averaged scales
avgSf = mean(Sf,[2 3]);
figure; plot(avgSf,'LineWidth',1.5);
xlim([0 255])
title('Spatial Average of $|Sf|$ with $\alpha = 0$ ','Interpreter','Latex','FontSize',14);
xlabel('$\tau$','Interpreter','Latex');
ylabel('Average $|Sf|$','Interpreter','Latex');
%figure; plot(avgSf/max(avgSf(:)));
hold on
qProceed = input('Proceed? (y/n) ');
if strcmp(qProceed,'n') 
    return
elseif strcmp(qProceed,'y')
    locs = input('Where to threshold? ');
if length(locs) == 1
    AllRegion = zeros(2,2);
    AllRegion(1,:) = [1 locs];
    AllRegion(2,:) = [locs+1 n];
else
    SizeRegion = length(locs);
    AllRegion = zeros(SizeRegion+1,2);
    AllRegion(1,:) = [1 locs(1)];
    for kk = 2:SizeRegion
        AllRegion(kk,:) = [locs(kk-1)+1, locs(kk)];
    end
    AllRegion(end,:) = [locs(end)+1 n];
end
for jj = 1:length(locs)
    xline(locs(jj))
end

return


% Matrix with all pixels to visit
% n^2 locations in order 
% [1 1], [1 2], ..., [1 256]
% [2 1, [2 2], ..., [2 256]
% etc 
[XX,YY] = meshgrid(1:256,1:256);
A=[XX(:),YY(:)];

[tM, tN] = size(f);

%tScaleMatrix = zeros(tM,tN,length(AllRegion));
tScaleMatrix = zeros(tM,tN);
row0 = 1;
col0 = 1;
% Check all pixels to see if t
% corresponding to max(Sf) is in Region N

[rowSizeRegions, ~] = size(AllRegion);

for ii = 1:length(XX)
    % update row
    rowCurrent = row0+(ii-1);
    for kk = 1:length(YY)
        colCurrent = col0+(kk-1);
        [~, indexOfMaxSf] = max(Sf(:,rowCurrent,colCurrent));
        tCurrent = tauAxis(indexOfMaxSf);
        % jj is the marker for which interval correspond to pixel 
%            for jj = 1:length(AllRegion)
        % Initialize jj
        for jj = 1:rowSizeRegions
           truthVal = tCurrent >= AllRegion(jj,1) && tCurrent <= AllRegion(jj,2);
           if truthVal == 1
               tScaleMatrix(rowCurrent,colCurrent) = jj;
           end
        end
    end
end

tRange = length(AllRegion);

figure;
imagesc(tScaleMatrix)
title('Region-wise Classification of Pixels','Interpreter','Latex');
axis off

%%%%%%%%%% End classification with maxima 

%%%%%%%%%% Begin classification with k-means 
SfRowVectors = zeros(M*N,n+1);

% Store "feature" row-vectors in a M*N (total pixels) by 256 (scales) matrix
for jj=1:M*N
    tempCnt = A(jj,:);
    tempSf = Sf(:,tempCnt(1),tempCnt(2));
%    temp = reshape(tempSf',1,theta_size*n);
    % Original
%    Sf_scales(jj,:) = Sf(:,:,tempCnt(1),tempCnt(2))';
    SfRowVectors(jj,:) = tempSf';
end

% AllRegion

% use cityblock metric for high dimensionality distance 
% look at other ways to cluster data: Hierarchical Clustering
% use PCA to reduce dimensionality 
[coeff,score,latent,tsquared,explained,mu] = pca(SfRowVectors);
varianceVec = []; kk = 1;
while sum(varianceVec(:)) < 95
    jj=explained(kk);
    varianceVec = [varianceVec jj];
    kk = kk+1;
end
PCA_end = kk;
reducedDimension = coeff;
reducedSfRowVectors = SfRowVectors*reducedDimension(:,1:PCA_end);

clstrs = 3;

[idx,C] = kmeans(reducedSfRowVectors,clstrs,'MaxIter',1000,'Replicate',3);
% Reshape index matrix into 256 by 256 
kmeans_index = reshape(idx, [M,N])';

figure;
imagesc(kmeans_index);
title('K-means Clustering','Interpreter','Latex');

max_clstrs = 10;

T = clusterdata(reducedSfRowVectors,'Linkage','ward','SaveMemory','on','Maxclust',max_clstrs);
agglo_cluster = reshape(T, [M,N])';
figure;
imagesc(agglo_cluster);
title('agglo clustering');
end

return 
%%
figure;
gryK = mat2gray(kmeans_index);
imshow(f);
green = cat(3, zeros(size(f)), ones(size(f)), zeros(size(f)));
hold on
h = imshow(green);
hold off
set(h,'AlphaData',gryK);

return 
%%
% read kmeans_cluster and agglo_cluster to proceed with the following
% figure
% they are matrices which store how pixels have been classified

figure;
subaxis(2,1,1, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0.1);
imagesc(k_idx);
axis off;
%title('\textit{k-Means} Clustering', 'Interpreter','Latex','FontSize',14);

subaxis(2,1,2, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0.1);
imagesc(aggglo_idx);
axis off;
%title('Hierarchical Clustering','Interpreter','Latex','FontSize',14);

C = colormap(cool(10));

return 
%%

%temp1 = load('AggoMatrix.mat'); 
temp2 = load('tScaleMatrix.mat');
temp3 = load('kMeansMatrix.mat');
temp4 = load('kMeansMatrix_3Clusters.mat');
%AggoMatrix = temp1.AggoMatrix;
tScaleMatrix = temp2.tScaleMatrix;
kMeansMatrix = temp3.kMeansMatrix;
kMeansMatrix_3Clusters = temp4.kMeansMatrix_3Clusters;

% Change figure 3.2 in Thesis to show original f and below it 
% the first type of classification 
figure;
f = yellowSquare(1);
ax1=subaxis(2,1,1, 'Spacing', 0.05, 'Padding', 0, 'Margin', 0.1);
%imshow(f);
f(:,1) = 0;
f(:,end) = 0;
f(1,:) = 0;
f(end,:) = 0;
imagesc(f);
colormap(ax1,gray)
title('Input Image $f$','Interpreter','Latex','FontSize',14);
axis off; 
ax2=subaxis(2,1,2, 'Spacing', 0.05, 'Padding', 0, 'Margin', 0.1);
imagesc(tScaleMatrix);
title('Classification Considering Regions $R_1$, $R_2$, and $R_3$','Interpreter','Latex','FontSize',14)
colormap(ax2,parula)
axis off;

sgtitle('Clustering of $|Sf|$ Based on Extrema in $R_i$','Interpreter','Latex','FontSize',14)
set(gcf, 'Color', 'w');
export_fig('3_RegionWiseClassification_Updated','-eps');

% Labels
% Play with obtained matrices to match color ... 

% Copy kMeansMatrix
kMean_Copy = kMeansMatrix;

% find indices which give a 6 and assign 10
Temp_Indices = find(kMeansMatrix == 6);
kMean_Copy(Temp_Indices) = 10;

% find indicies which give a 10 and assign 6 
Temp_Indices = find(kMeansMatrix == 10);
kMean_Copy(Temp_Indices) = 6;

% find indices which give a 4 and assign 1 
Temp_Indices = find(kMeansMatrix == 4);
kMean_Copy(Temp_Indices) = 1;

% find indices which give a 1 and assign 4
Temp_Indices = find(kMeansMatrix == 1);
kMean_Copy(Temp_Indices) = 4;

% region-wise clustering 
figure;
subaxis(2,1,1, 'Spacing', 0.05, 'Padding', 0, 'Margin', 0.1);
imagesc(tScaleMatrix);
%axis off;
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
title('Region-wise Classification','Interpreter','Latex')
set(gcf, 'Color', 'w');

% kmeans clustering
subaxis(2,1,2, 'Spacing', 0.05, 'Padding', 0, 'Margin', 0.1);
imagesc(kMean_Copy);
title('Clustering via \textit{k-means}','Interpreter','Latex')
%axis off
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
set(gcf, 'Color', 'w');

sgtitle('Comparison of Clustering Techniques','Interpreter','Latex','FontSize',14)
export_fig('3_kMeansAndHierarchichal','-eps');

Indices = find(kMeansMatrix == 10)
kMean_Copy(Indices) = K 

%% recreate graphs
% R_i Method vs kMeans (3 clusters)
temp2 = load('tScaleMatrix.mat');
temp3 = load('kMeansMatrix.mat');
%AggoMatrix = temp1.AggoMatrix;
tScaleMatrix = temp2.tScaleMatrix;
kMeansMatrix = temp3.kMeansMatrix;

Ri_Mat = tScaleMatrix;

figure;
subaxis(2,1,1, 'Spacing', 0.05, 'Padding', 0, 'Margin', 0.1);
imagesc(Ri_Mat);
%axis off;
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
title('$R_i$ Method','Interpreter','Latex')
set(gcf, 'Color', 'w');

% kmeans clustering
kMeans_3Clusters = kmeans_3;
subaxis(2,1,2, 'Spacing', 0.05, 'Padding', 0, 'Margin', 0.1);
imagesc(kMeans_3Clusters);
title('Three Clusters via \textit{k-means}','Interpreter','Latex')
%axis off
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
set(gcf, 'Color', 'w');

sgtitle('Comparison of Clustering Techniques','Interpreter','Latex','FontSize',14)
export_fig('3_RiMethodVSkMeans','-eps');

%% recreate graphs 
% kMeans (3 clusters) vs kMeans (10 clusters)
temp2 = load('tScaleMatrix.mat');
temp3 = load('kMeansMatrix.mat');
tScaleMatrix = temp2.tScaleMatrix;
kMeansMatrix = temp3.kMeansMatrix;

kMeans_3Clusters = kmeans_3;
figure;
subaxis(2,1,1, 'Spacing', 0.05, 'Padding', 0, 'Margin', 0.1);
imagesc(kMeans_3Clusters);
%axis off;
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
title('Three Clusters','Interpreter','Latex')
set(gcf, 'Color', 'w');

% kmeans clustering
%kMeans_10Clusters = kMeansMatrix;
kMeans_10Clusters = kMean_Copy;
subaxis(2,1,2, 'Spacing', 0.05, 'Padding', 0, 'Margin', 0.1);
imagesc(kMeans_10Clusters);
title('Ten Clusters','Interpreter','Latex')
%axis off
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
set(gcf, 'Color', 'w');

sgtitle('Comparison of \textit{k-means}','Interpreter','Latex','FontSize',14)
export_fig('3_kMeansComparison','-eps');

%%
kmeans_3 = 0*kmeans_index;
% find indices which give a 2 and assign 3
Temp_Indices = find(kmeans_index == 2);
kmeans_3(Temp_Indices) = 3;

Temp_Indices = find(kmeans_index == 3);
kmeans_3(Temp_Indices) = 2;

Temp_Indices = find(kmeans_index == 1);
kmeans_3(Temp_Indices) = 1;