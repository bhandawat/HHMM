function [closestPtBef,closestPtAft,idxB] = findBeforeAfter(A,B,type)
% This function finds the closest value of B to each value of A. 
%
% Inputs:
%    A: A matrix of values
%    B: A matrix of reference values
%    Type: Type of points we want to return (3 options)
%       'both' = find the closest point of B that is smaller than A and the
%       closest point that is larger than A.
%       'before' = find the closest point of B that is smaller than A.
%       'after' = find the closest point of B that is larger than A.

% Outputs:
%    closestPtBef: The values of B that is closest to and smaller than each
%           point of A (same dimensions as A)
%    closestPtAft: The values of B that is closest to and larger than each
%           point of A (same dimensions as A)
%    idxB: The index values of B for closestPtBef and closestPtAfter

[m,n] = size(A);
[B,ndxB] = sort(reshape(B,[],1));
A = reshape(A,[],1);

TMP = bsxfun(@(x,y) abs(x-y), A(:), B');
[~, idxB] = min(TMP,[],2);
closestPt = B(idxB);


if strcmpi(type,'both')
    direction = closestPt-reshape(A,[],1);

    idxBBef = idxB;
    idxBAft = idxB;

    idxBBef(direction>-1) = idxBBef(direction>-1)-1;
    idxBAft(direction<1) = idxBAft(direction<1)+1;

    BadPtsBef = idxBBef==0 | idxBBef>numel(B);
    BadPtsAft = idxBAft==0 | idxBAft>numel(B);
    
    idxBBef(BadPtsBef) = 1;idxBAft(BadPtsAft) = 1;
    
    closestPtBef = B(idxBBef)';
    closestPtAft = B(idxBAft)';
    closestPtBef(BadPtsBef) = NaN;
    closestPtAft(BadPtsAft) = NaN;
    
    closestPtBef = reshape(closestPtBef,m,n);
    closestPtAft = reshape(closestPtAft,m,n);
    
    idxB = [ndxB(idxBBef),ndxB(idxBAft)];
    
elseif strcmpi(type,'before')
    idxBBef = idxB;
    direction = closestPt-A;
    idxBBef(direction>-1) = idxBBef(direction>-1)-1;
    
    BadPts = idxBBef==0 | idxBBef>length(B);
    idxBBef(BadPts) = 1;
    closestPtBef = B(idxBBef)';
    closestPtBef(BadPts) = NaN;
    closestPtBef = reshape(closestPtBef,m,n);
    closestPtAft = [];
    
    idxB = ndxB(idxBBef);
    
elseif strcmpi(type,'after')
    idxBAft = idxB;
    direction = closestPt-A;
    idxBAft(direction<1) = idxBAft(direction<1)+1;
    
    BadPts = idxBAft==0 | idxBAft>length(B);
    idxBAft(BadPts) = 1;
    closestPtAft = B(idxBAft)';
    closestPtAft(BadPts) = NaN;
    closestPtAft = reshape(closestPtAft,m,n);
    closestPtBef = [];
    
    idxB = ndxB(idxBAft);
    
end

end