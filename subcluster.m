function [centers,sigmas,maxPoint,X,potVals] = subcluster(X,RA,xBounds,options)
% SUBCLUSTER Locates data cluster centers using subtractive clustering.
%
%   SUBCLUSTER finds the optimal data point to define a cluster
%   center based on the density of surrounding data points. All data points
%   within the distance RA of this point are then removed, in order to
%   determine the next data cluster and its center. This process is repeated
%   until all of the data is within the distance RA of a cluster center.
%
%   [C] = SUBCLUSTER(X,RA) clusters the data points in the M-by-N matrix X,
%   where M is the number of data points and N is the number of coordinates
%   (data dimensions) that specify each point. RA has a value between 0 and 1
%   and specifies the size of the cluster in each of the data dimensions,
%   assuming the data fall within a unit hyperbox (range [0 1]). Specifying a
%   smaller cluster radius will usually yield more (smaller) clusters in the
%   data. When RA is a scalar it is applied to all data dimensions, when it
%   is a vector, it has one entry for each data dimension. The cluster centers
%   are returned as rows of the matrix C. C has size J-by-N, if J clusters are
%   required to cluster the data.
%
%   [C] = SUBCLUSTER(...,XBOUNDS) also specifies a matrix XBOUNDS, of size
%   2-by-N, used to normalize the data X into a unit hyperbox (range [0 1]).
%   Each column of XBOUNDS provides the minimum and maximum values for the
%   corresponding data dimension. If XBOUNDS is an empty matrix or not provided,
%   the minimum and maximum data values found in X, are used as defaults.
%
%   [C] = SUBCLUSTER(...,OPTIONS) specifies a vector for changing the default
%   algorithm parameters:
%      OPTIONS(1):  The squash factor, is used to multiply the RA values to
%                   determine the neighborhood of a cluster center within which
%                   the existence of other cluster centers are discouraged.
%      OPTIONS(2):  The accept ratio, sets the potential, as a fraction of the
%                   potential of the first cluster center, above which another
%                   data point will be accepted as a cluster center.
%      OPTIONS(3):  The reject ratio sets the potential, as a fraction of the
%                   potential of the first cluster center, below which a data
%                   point will be rejected as a cluster center.
%      OPTIONS(4):  Displays progress information unless it is set to zero.
%   The default values for the OPTIONS vector are [1.25 0.5 0.15 0].
%
%   [C,S] = SUBCLUSTER(X,RA,XBOUNDS,OPTIONS) returns the cluster centers in the matrix C;
%   each row of C contains the position of a cluster center.
%   The returned S vector contains the sigma values that specify
%   the range of influence of a cluster center in each of the data dimensions.
%   All cluster centers share the same set of sigma values.
%
%   Examples
%       X1 = 10*rand(600,1);
%       X2 =  5*rand(600,1);
%       X3 = 20*rand(600,1)-10;
%       X = [X1 X2 X3];
%      [C] = subcluster(X,0.5) specifies a range of influence of 0.5 for all data
%      dimensions.
%
%      [C] = subcluster(X,[0.5 0.25 0.3],[],[2.0 0.8 0.7 0]) specifies a range of
%      influence of 0.5, 0.25, and 0.3 for the first, second and third data
%      dimensions. The scaling factors for mapping the data into a unit hyperbox
%      will be obtained from the minimum and maximum data values. The squash
%      factor is set to 2.0, indicating that we want to only find clusters that
%      are far from each other, the accept ratio is set to 0.8, indicating that
%      we will only accept data points that have very strong potential of being
%      cluster centers, the reject ratio is set to 0.7, indicating that we want
%      to reject all data points without a strong potential.
%

%	Reference 
%	S. Chiu, "Fuzzy Model Identification Based on Cluster Estimation," J. of
%	Intelligent & Fuzzy Systems, Vol. 2, No. 3, 1994.

[q,n] = size(X);
if nargin < 4
	options = [1.25 0.5 0.15 0]; % [squash factor,accept ratio,reject ratio,displays progress]
end

if nargin < 3
    xBounds = [];
end

% if only one value is given as the range of influence, then apply that
% value to all data dimensions 
if length(RA) == 1 & n ~= 1
    RA = RA * ones(1,n);
end

sqshFactor  = options(1);
acceptRatio = options(2);
rejectRatio = options(3);
verbose     = options(4);

if verbose
    disp('Normalizing data...');
end

if length(xBounds) == 0
    % No data scaling range values are specified, use the actual minimum and maximum values
    minX = min(X);
    maxX = max(X);
    % If the actual min and max values have a range of zero, calculate a small range
    % relative to the data, to allow a sigma value to be calculated.
    index = find(maxX == minX);
    minX(index) = minX(index) - 0.0001*(1 + abs(minX(index)));
    maxX(index) = maxX(index) + 0.0001*(1 + abs(maxX(index)));
else
    % Use the user supplied data range values in xBounds
    minX = xBounds(1,:);
    maxX = xBounds(2,:);
    % Verify correct dimensions and values for xBounds were supplied
    if length(minX) ~=  size(X,2)
        error('xBounds contains the wrong dimensions for input data X');
    elseif any(maxX == minX)
        error('xBounds has a data range of zero');
    end
end

% Normalize the data into a unit hyperbox using the verified minX and maxX
for p=1:q,
    for i = 1:n
        X(p,i) = (X(p,i) - minX(i)) / (maxX(i) - minX(i));
    end
end
% X = min(max(X,0),1);

if verbose
    disp('Computing potential for each data point...');
end

% potVals = the potential of each data point to be a cluster center
% potVals = zeros(1,q);

% compute the initial potentials for each data point

UseWaitBar = (q > 500);
if UseWaitBar
	h = waitbar(0,'Computing clusters. Please wait...');
end

for p=1:q,
    potVals(p) = 0.0;
    for pp=1:q,
        d2(pp) = 0.0;
        for i=1:n,
            dx2(pp,i) = (X(p,i)-X(pp,i))/RA(i);
            d2(pp) = d2(pp)+dx2(pp,i)^2;
        end
        potVals(p) = potVals(p) + exp(-4*d2(pp));
    end
        
	if UseWaitBar & mod(p,100)==0
		waitbar(p/q);
    end
end

if UseWaitBar
	close(h)    
end

% Find the data point with highest potential value.  refPotVal is the
% highest potential value, used as a reference for accepting/rejecting
% other data points as cluster centers.
[refPotVal,maxPotIndex] = max(potVals);

% Start iteratively finding cluster centers and subtracting potential
% from neighboring data points.  maxPotVal is the current highest
% potential value and maxPotIndex is the associated data point's index.
maxPotVal = refPotVal;

% centers = the cluster centers that has been found
numClusters = 0;
findMore = 1;

while findMore & maxPotVal
	findMore = 0;
    for i=1:n,
        maxPoint(i) = X(maxPotIndex,i);
    end
	maxPotRatio = maxPotVal/refPotVal;
	if maxPotRatio > acceptRatio
		% the new peak value is significant, accept it
		findMore = 1;
	elseif maxPotRatio > rejectRatio
		% accept this data point only if it achieves a good balance between having
		% a reasonable potential and being far from all existing cluster centers
		minDistSq = -1;
		for j=1:numClusters
            dxSq = 0.0;
            for i=1:n,
                dx1(i) = (maxPoint(i) - centers(j,i))/RA(i);
                dxSq = dxSq + dx1(i)^2;
            end
			if minDistSq < 0 | dxSq < minDistSq
				minDistSq = dxSq;
            end
        end
		minDist = sqrt(minDistSq);
		if (maxPotRatio + minDist) >= 1
			findMore = 1;	% tentatively accept this data point as a cluster center
		else
			findMore = 2;	% remove this point from further consideration, and continue
		end
	end	% end if maxPotRatio > acceptRatio
	if findMore == 1
		% add the data point to the list of cluster centers
		numClusters = numClusters + 1;
        for i=1:n,
            centers(numClusters,i) = maxPoint(i);
        end
		if verbose
			msg = sprintf('Found cluster %g, potential = %g',numClusters,maxPotRatio);
			disp(msg);
        end
		% subtract potential from data points near the new cluster center
		for pp=1:q,
            d2 = 0.0;
            for i=1:n,
                dx2(pp,i) = (maxPoint(i)-X(pp,i))/(sqshFactor*RA(i));
                d2 = d2 + dx2(pp,i)^2;
            end
            deduct(pp) = maxPotVal*exp(-4*d2);            
            potVals(pp) = potVals(pp) - deduct(pp);
            if potVals(pp) < 0,
                potVals(pp) = 0;
            end
        end        
		% find the data point with the highest remaining potential
		[maxPotVal,maxPotIndex] = max(potVals);
	elseif findMore == 2
		potVals(maxPotIndex) = 0;
		[maxPotVal,maxPotIndex] = max(potVals);
	end % end if findMore == 1
end % end while findMore & maxPotVal

% Scale the cluster centers from the normalized values back to values in
% the original range
for j=1:numClusters,
    for i=1:n,
        centers(j,i) = (centers(j,i) * (maxX(i) - minX(i))) + minX(i);
    end
end
% Compute the sigma values for the clusters
for i=1:n,
    sigmas(i) = (RA(i) .* (maxX(i) - minX(i))) / sqrt(8.0);
end
