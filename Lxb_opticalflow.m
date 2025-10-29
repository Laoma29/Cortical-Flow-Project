function [flowField, int_dF, errorData, errorReg, poincare, interval] = ...
  bst_opticalflow(F, FV, Time, tStart, tEnd, hornSchunck,sample_rate)
% BST_OPTICALFLOW: Computes the optical flow of MEG/EEG activities on the cortical surface.
%% USAGE:  [flow, dEnergy, int_dI, errorData, errorReg, index] = ...
%   bst_opticalflow(dataFile, channelFile, tStart, tEnd, hornSchunck)
%
% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2020 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Julien Lef�vre, 2006-2010
%          Syed Ashrafulla, 2010
％　　　　Ｘｉａｏｂｏ　Ｌｉｕ　２０２０

dimension = 3; % 2 for projected maps
Faces = FV.Faces; Vertices = FV.Vertices; VertNormals = FV.VertNormals;
nVertices = size(Vertices,1); % VertNormals = FV.VertNormals';
nFaces = size(Faces,1);
tStartIndex = find(Time < tStart-100*eps, 1, 'last')+1; % Index of first time point for flow calculation
if isempty(tStartIndex)
    [tmp, tStartIndex] = min(Time);
    tStartIndex = tStartIndex + 1;
end
tEndIndex = find(Time < tEnd-100*eps, 1, 'last')+1; % Index of last time point for flow calculation
if isempty(tEndIndex)
    [tmp, tEndIndex] = max(Time);
    tEndIndex = tEndIndex + 1;
end
if tEndIndex > size(Time, 2)
    tEndIndex = size(Time, 2);
end
interval = Time(tStartIndex:tEndIndex); % Interval of flow calculations
intervalLength = tEndIndex-tStartIndex+1; % Length of time interval for calculations
M = max(max(abs(F(:,tStartIndex-1:tEndIndex)))); F = F/M; % Scale values to avoid precision error
flowField = zeros(nVertices, dimension, intervalLength);
dEnergy = zeros(1, intervalLength);
errorData = zeros(1, intervalLength);
errorReg = zeros(1, intervalLength);
int_dF = zeros(1, intervalLength);

[regularizerOld, gradientBasis, tangentPlaneBasisCell, tangentPlaneBasis, ...
  triangleAreas, FaceNormals, sparseIndex1, sparseIndex2]= ...
  regularizing_matrix(Faces, Vertices, VertNormals, dimension); % 2 for projected maps
regularizer = hornSchunck * regularizerOld;

% Projection of flow on triangle (for Poincar� index)
Pn = zeros(3,3,nFaces);
parfor facesIdx = 1:nFaces
    Pn(:,:,facesIdx) = eye(3) - ...
      (FaceNormals(facesIdx,:)'*FaceNormals(facesIdx,:));
end
poincare = zeros(nFaces, intervalLength);

% Optical flow calculations
bst_progress('start', 'Optical Flow', ...
  'Computing optical flow ...', tStartIndex-1, tEndIndex);
for timePoint = tStartIndex:tEndIndex
  timeIdx = timePoint-tStartIndex+1;
  
  % Solve for flow gpuArray
  dF = (F(:,timePoint)-F(:,timePoint-1))/(1/sample_rate);%;
%   dF = F(:,timePoint)-F(:,timePoint-1);

  [dataFit, B] = data_fit_matrix(Faces, nVertices, ...
    gradientBasis, triangleAreas, tangentPlaneBasisCell, ...
    F(:,timePoint), dF, dimension, ...
    sparseIndex1, sparseIndex2);
  X = (dataFit + regularizer) \ B;
  
  % Save flows correctly (for 3D surface, have to project back to 3D)
  if dimension == 3 % X coordinates are in tangent space
      flowField(:,:,timeIdx) = ...
        repmat(X(1:end/2), [1,3]) .* tangentPlaneBasis(:,:,1) ...
        + repmat(X(end/2+1:end), [1,3]) .* tangentPlaneBasis(:,:,2);
  else % X coordinates are in R^2
      flowField(:,:,timeIdx) = [X zeros(nVertices,1)];
  end
  
  errorData(timeIdx)= X'*dataFit*X - 2*B'*X; % Data fit error
  errorReg(timeIdx)= X'*regularizer*X; % Regularization error

  % Variational formulation constant
  dF12=(dF(Faces(:,1),:)+dF(Faces(:,2))).^2;
  dF23=(dF(Faces(:,2),:)+dF(Faces(:,3))).^2;
  dF13=(dF(Faces(:,1),:)+dF(Faces(:,3))).^2;
  int_dF(timeIdx) = sum(triangleAreas.*(dF12+dF23+dF13)) / 24;
  
  % Precompute flowfield with faces to save time in the loop.
  FacesFlowField = reshape(flowField(Faces', :, timeIdx)', [3,3,nFaces]);
  
  % Poincare Index of each triangle
  parfor facesIdx=1:nFaces
      poincare(facesIdx,timeIdx) = ... % projection of flowField(f,:,t) on triangle f
          poincare_index(Pn(:,:,facesIdx) * FacesFlowField(:,:,facesIdx));
  end

  % Displacement energy
  v12 = sum((flowField(Faces(:,1),:,timeIdx)+flowField(Faces(:,2),:,timeIdx)).^2,2) / 4;
  v23 = sum((flowField(Faces(:,2),:,timeIdx)+flowField(Faces(:,3),:,timeIdx)).^2,2) / 4;
  v13 = sum((flowField(Faces(:,1),:,timeIdx)+flowField(Faces(:,3),:,timeIdx)).^2,2) / 4;
  dEnergy(timeIdx) = sum(triangleAreas.*(v12+v23+v13));
  
  bst_progress('inc', 1); % Update progress bar
end 

bst_progress('stop');
  
end
