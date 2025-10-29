
%% ===== TESSELATION TANGENT BUNDLE =====
function tangentPlaneBasis = basis_vertices(VertNormals)

nVertices = size(VertNormals,1); VertNormals = -VertNormals;
tangentPlaneBasis = zeros(nVertices, 3, 2); % Initialize

tangentPlaneBasis(:,:,1) = diff(VertNormals(:, [2 3 1 2]).').';

bad = abs(dot(VertNormals, ones(nVertices,3)/sqrt(3), 2)) > 0.97;
tangentPlaneBasis(bad,:,1) = ...
    [VertNormals(bad,2) -VertNormals(bad,1) zeros(sum(bad),1)];

tangentPlaneBasis(:,:,2) = cross(VertNormals, tangentPlaneBasis(:,:,1));


tangentPlaneBasisNorm = sqrt(sum(tangentPlaneBasis.^2,2));
tangentPlaneBasis = tangentPlaneBasis ./ tangentPlaneBasisNorm(:,[1 1 1],:);

end
