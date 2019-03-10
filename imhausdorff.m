function [D,idx]=imhausdorff(A,B,method)
%IMHAUSDORFF distance for image segmentation (i.e. of gridded data).
%
%function [D,IDX]=imhausdorff(A,B,method='euclidean')
%
%   D = IMHAUSDORFF(BW1,BW2) computes the Hausdorff distance of binary
%   images BW1 and BW2
%
%   D = IMHAUSDORFF(L1,L2) computes the Hausdorff distance for each label
%   in label images L1 and L2.
%
%   D = IMHAUSDORFF(C1,C2) computes the Hausdorff distance for each
%   category in categorical images C1 and C2.
%
%   D = IMHAUSDORFF(A,B,METHOD) lets you compute the Hausdorff distance
%   with an alternate point-to-point distance.  METHOD can be
%   'cityblock', 'chessboard', 'quasi-euclidean', 'euclidean' or
%   'euclidean_precise'.  METHOD defaults to 'euclidean' if not specified.
%   See METHOD in BWDIST.
%   'euclidean_precise' aims at double precision by re-evaluating the
%   distance of the farthest points found using the closest-pixel map of
%   BWDIST
%
%   [D,IDX] = IMHAUSDORFF(...) also returns an N x 2 matrix of linear
%   indices of the resp. farthest points contributing to the distance.
%
%   Notes
%   -----
%   IMHAUSDORFF uses BWDIST for fast computation, therefore the output
%   is only of single precision unless METHOD is 'euclidean_precise'.
%
%   IMHAUSDORFF(A,[]) = inf, nnz(A)>0
%   IMHAUSDORFF([],[]) = 0
%
%   Example 1
%   ---------
%   Compute the Hausdorff distance for a binary segmentation.
%
%     % Read in an image with an object we wish to segment.
%     A = imread('hands1.jpg');
%
%     % Convert the image to grayscale.
%     I = rgb2gray(A);
%
%     % Use active contours to segment the hand.
%     mask = false(size(I));
%     mask(25:end-25,25:end-25) = true;
%     BW = activecontour(I, mask, 300);
%
%     % Read in the ground truth against which to compare the segmentation.
%     BW_groundTruth = imread('hands1-mask.png');
%
%     % Compute the Hausdorff distance of this segmentation.
%     [distance,idx] = imhausdorff(BW, BW_groundTruth);
%
%     % Display both masks on top of one another.
%     figure
%     imshowpair(BW, BW_groundTruth)
%     title(['Hausdorff distance = ' num2str(distance)])
%
%     % Display a line indicating the farthest points
%     [y,x]=ind2sub(size(BW),idx);
%     hold on;
%     plot(x,y,'rx-','linewidth',2');
%     hold off
%
%   Example 2
%   ---------
%   Compute the Hausdorff distance for for n-ary segmentation.
%
%     % Read in an image with several objects we wish to segment.
%     RGB = imread('yellowlily.jpg');
%
%     % Create scribbles for three regions.
%     region1 = [350 700 425 120]; % [x y w h] format
%     BW1 = false(size(RGB,1),size(RGB,2));
%     BW1(region1(2):region1(2)+region1(4),region1(1):region1(1)+region1(3)) = true;
%  
%     region2 = [800 1124 120 230];
%     BW2 = false(size(RGB,1),size(RGB,2));
%     BW2(region2(2):region2(2)+region2(4),region2(1):region2(1)+region2(3)) = true;
%  
%     region3 = [20 1320 480 200; 1010 290 180 240]; 
%     BW3 = false(size(RGB,1),size(RGB,2));
%     BW3(region3(1,2):region3(1,2)+region3(1,4),region3(1,1):region3(1,1)+region3(1,3)) = true;
%     BW3(region3(2,2):region3(2,2)+region3(2,4),region3(2,1):region3(2,1)+region3(2,3)) = true;
%
%     % Display the seed regions on top of the image.
%     figure
%     imshow(RGB)
%     hold on
%     visboundaries(BW1,'Color','r');
%     visboundaries(BW2,'Color','g');
%     visboundaries(BW3,'Color','b');
%     title('Seed regions')
%
%     % Segment the image into three regions using
%     % geodesic distance-based color segmentation.
%     L = imseggeodesic(RGB,BW1,BW2,BW3,'AdaptiveChannelWeighting',true);
%
%     % Load a ground truth segmentation of the image into three regions.
%     L_groundTruth = double(imread('yellowlily-segmented.png'));
%
%     % Visually compare the segmentation results with the ground truth.
%     figure
%     imshowpair(label2rgb(L),label2rgb(L_groundTruth),'montage')
%     title('Comparison of segmentation results (left) and ground truth (right)')
%
%     % Compute the Hausdorff distance
%     % for each segmented region.
%     [distance,idx] = imhausdorff(L, L_groundTruth);
%     distance
%
%     % Display a line indicating the farthest points
%     [y,x]=ind2sub(size(L),idx.');
%     hold on;
%     plot(x,y,'kx-','linewidth',2');
%     text(mean(x),mean(y),num2str(distance),'VerticalAlignment','bottom');
%     hold off
%
%   See also HAUSDORFF, BWDIST, BFSCORE, DICE, JACCARD.
%
%Author: Joakim Lindblad

% Copyright (c) 2019, Joakim Lindblad


validateInput = @(x,name,pos) validateattributes(x, ...
    {'logical','numeric','categorical'}, ...
    {'real','nonempty','nonsparse'}, ...
    mfilename,name,pos);

validateInput(A,'A',1);
validateInput(B,'B',2);

if any(size(A) ~= size(B))
    error(message('images:validate:unequalSizeMatrices','A','B'))
end

if ~isa(A,class(B))
    error(message('images:validate:differentClassMatrices','A','B'))
end

if isa(A,'numeric')
    % Additional validation for label matrices
    validateLabelMatrix = @(x,name,pos) validateattributes(x, ...
        {'numeric'}, ...
        {'finite','nonnegative','integer'}, ...
        mfilename,name,pos);
    validateLabelMatrix(A,'L1',1);
    validateLabelMatrix(B,'L2',2);
end

if isa(A,'categorical')
    % Additional validation for categorical matrices
    if ~isequal(categories(A),categories(B))
        error(message('images:segmentation:nonIdenticalCategories2', ...
            'C1','C2'))
    end
end

if nargin < 3
    method = 'euclidean';
else
    valid_methods = {'euclidean','euclidean_precise','cityblock', ...
        'chessboard','chebychev','quasi-euclidean'};
    method = validatestring(method, valid_methods, ...
                          mfilename, 'METHOD', 2);
    if strcmp(method,'chebychev')
        method = 'chessboard'; % synonymous
    end
end


if isa(A,'logical')
    % binary segmentation
    if (nargout == 2)
        [D,idx] = bwhdist2(A,B,method);
    else
        D = bwhdist(A,B,method);
    end
else
    if isa(A,'categorical')
        % categorical matrices
        classes = categories(A);
    else
        % label matrices
        classes = (1:max(max(A(:)),max(B(:))))';
    end
    
    BW1 = convertToCellOfLogicals(A,classes);
    BW2 = convertToCellOfLogicals(B,classes);

    if strcmp(method,'euclidean_precise')
        D = nan(numel(classes),1,'double');
    else
        D = nan(numel(classes),1,'single');
    end
    if (nargout == 2)
        idx=nan(numel(classes),2);
        for k = 1:numel(classes)
            [D(k),idx(k,:)] = bwhdist2(BW1{k},BW2{k},method);
        end
    else
        for k = 1:numel(classes)
            D(k) = bwhdist(BW1{k},BW2{k},method);
        end
    end
end


% Max of distance from A to B and from B to A
function D=bwhdist(A,B,method)
D=max( dir_hdist(A,B,method), dir_hdist(B,A,method) );

% Max distance from points in A to B
function D=dir_hdist(A,B,method)
if strcmp(method,'euclidean_precise')
    D=dir_hdist_euclidean_precise(A,B);
    return;
end
if any(A(:))
    DT=bwdist(B,method);
    D=max(DT(A));
else
    D=single(0); % Distance from [] to B is always zero
end


% Max of distance from A to B and from B to A
% Also returning linear indices of fathest point pair
function [D,idx]=bwhdist2(A,B,method)
[d_AB,idx_AB]=dir_hdist2(A,B,method);
[d_BA,idx_BA]=dir_hdist2(B,A,method);
if(d_AB>d_BA)
    D=d_AB;
    idx=idx_AB;
else
    D=d_BA;
    idx=fliplr(idx_BA); % First idx is in A
end

% Max distance from points in A to B
% Also returning linear indices of fathest point pair
function [D,idx]=dir_hdist2(A,B,method)
if strcmp(method,'euclidean_precise')
    [D,idx]=dir_hdist_euclidean_precise(A,B);
    return;
end
if any(A(:))
    [DT,idxB]=bwdist(B,method);
    A=find(A); % Set points in A
    idxB=idxB(A); % Only care about points under A
    DT=DT(A);
    [D,idxA]=max(DT); % Farthest point in A (out of A~=0)
    idxB=idxB(idxA); % Used point in B

    if idxB>0
        idxA=A(idxA); % Used point in A
        idx=[idxA,idxB];
    else
        idx=[0,0];
    end
else
    D=single(0); % Distance from [] to B is always zero
    idx=[0,0];
end


% Precise max Euclidean distance from every point in A to B
% Re-computing distance for the found maximum distance pair
function [D,idx]=dir_hdist_euclidean_precise(A,B)
if any(A(:))
    [DT,idxB]=bwdist(B,'euclidean');
    A=find(A); % Set points in A
    idxB=idxB(A); % Only care about points under A
    DT=DT(A);
    [~,idxA]=max(DT); % Farthest point in A (out of A~=0)
    idxB=idxB(idxA); % Used point in B

    if idxB>0
        idxA=A(idxA); % Used point in A
        idx=[idxA,idxB];
        s=cell(1,ndims(B));
        [s{:}]=ind2sub(size(B),[idxA;idxB]); % Subscripts of the used points
        D=norm(diff([s{:}],[],1),2); % Euclidean distance
    else
        D=inf;
        idx=[0,0];
    end
else
    D=0; % Distance from [] to B is always zero
    idx=[0,0];
end


% One binary image for each label
function C=convertToCellOfLogicals(A,classes)
C=arrayfun(@(k) A==k, classes ,'UniformOutput', false);
