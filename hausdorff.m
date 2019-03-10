function HD = hausdorff(a,b,method)
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
%   See also IMHAUSDORFF, PDIST2.
%
%Author: Joakim Lindblad

% Copyright (c) 2019, Joakim Lindblad

if isempty(a) || isempty(b)
	if isempty(a) && isempty(b)
		HD=0;
	else
		HD=inf;
	end
	return
end

% Max of dist from A to B and B to A
if (size(a,1)*size(b,1)) < 256*256
	D = pdist2(a,b,'squaredeuclidean');
	HD = sqrt(max(max(min(D,[],1)),max(min(D,[],2))));
else
	% Less memory hungry version
	D1 = pdist2(a,b,'squaredeuclidean','Smallest',1);
	D2 = pdist2(b,a,'squaredeuclidean','Smallest',1);
	HD = sqrt(max(max(D1),max(D2)));
end
