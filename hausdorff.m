function [D,idx] = hausdorff(A,B,method)
%HAUSDORFF distance for point clouds.
%
%function [D,IDX]=hausdorff(a,b,method='euclidean')
%
%   D = HAUSDORFF(A,B) computes the Hausdorff distance between
%   point sets A and B. Rows of A and B correspond to observations,
%   and columns correspond to variables. A and B must have same number
%   of columns.
%
%   D = HAUSDORFF(A,B,METHOD) lets you compute the Hausdorff distance
%   with an alternate point-to-point distance.  METHOD can be any
%   method supported by PDIST2.  METHOD defaults to 'euclidean' if not
%   specified.
%
%   D = HAUSDORFF(A,B,DISTFUN) lets you compute the Hausdorff distance
%   with a distance function specified using a function handle @
%
%   [D,IDX] = HAUSDORFF(...) also returns the indices of the farthest
%   points contributing to the distance.
%
%   Notes
%   -----
%   HAUSDORFF uses PDIST2 for computation. For gridded image data
%   it is often preferred to use IMHAUSDORFF.
%
%   HAUSDORFF(A,[]) = inf
%   HAUSDORFF([],[]) = 0
%
%   Example 1
%   ---------
%   Compute the Hausdorff distance for a binary segmentation.
%   This is much slower than using IMHAUSDORFF.
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
%     % Extract object point coordinates
%     A=regionprops(BW,'PixelList');
%     B=regionprops(BW_groundTruth,'PixelList');
%
%     % Compute the Hausdorff distance of this segmentation.
%     [distance,idx] = hausdorff(A.PixelList,B.PixelList);
%
%     % Display both masks on top of one another.
%     figure
%     imshowpair(BW, BW_groundTruth)
%     title(['Hausdorff distance = ' num2str(distance)])
%
%     % Display a line indicating the farthest points
%     p=[A.PixelList(idx(1),:);B.PixelList(idx(2),:)];
%     hold on;
%     plot(p(:,1),p(:,2),'rx-','linewidth',2');
%     hold off
%
%   See also IMHAUSDORFF, PDIST2.
%
%Author: Joakim Lindblad

% Copyright (c) 2019, Joakim Lindblad


if size(A,2) ~= size(B,2)
    error('A and B must have the same number of columns.');
end

if nargin < 3
    method = 'euclidean';
else
    if strcmp(method,'chessboard')
        method = 'chebychev'; % synonymous
    end
end


if isempty(A) || isempty(B)
    if isempty(A) && isempty(B)
        HD=0;
    else
        HD=inf;
    end
    return
end

if strcmp(method,'euclidean')
    method='squaredeuclidean'; % faster
    apply_root=true;
else
    apply_root=false;
end


% Max of dist from A to B and B to A
if (size(A,1)*size(B,1) < 1e8)
    D = pdist2(A,B,method);
    [D1,idxA1] = min(D,[],1);
    [D2,idxB1] = min(D,[],2);
    clear D;
else
    % Less memory hungry version
    [D1,idxA1] = pdist2(A,B,method,'Smallest',1);
    [D2,idxB1] = pdist2(B,A,method,'Smallest',1);
end
[D1,idxB2]=max(D1);
[D2,idxA2]=max(D2);

if (D1>D2)
    D=D1;
    idx=[idxA1(idxB2),idxB2];
else
    D=D2;
    idx=[idxA2,idxB1(idxA2)];
end

if apply_root
    D = sqrt(D);
end

if (nargout < 2)
    clear idx;
end
