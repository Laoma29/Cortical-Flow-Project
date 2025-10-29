
function [gradientBasis, triangleAreas, FaceNormals] = ...
  geometry_tesselation(Faces, Vertices, dimension)

u = Vertices(Faces(:,2),:)-Vertices(Faces(:,1),:);
v = Vertices(Faces(:,3),:)-Vertices(Faces(:,2),:);
w = Vertices(Faces(:,1),:)-Vertices(Faces(:,3),:);


uu = sum(u.^2,2);
vv = sum(v.^2,2);
ww = sum(w.^2,2);
uv = sum(u.*v,2);
vw = sum(v.*w,2);
wu = sum(w.*u,2);


h1 = w-((vw./vv)*ones(1,dimension)).*v;
h2 = u-((wu./ww)*ones(1,dimension)).*w;
h3 = v-((uv./uu)*ones(1,dimension)).*u;
hh1 = sum(h1.^2,2);
hh2 = sum(h2.^2,2);
hh3 = sum(h3.^2,2);


gradientBasis = cell(1,dimension);
gradientBasis{1} = h1./(hh1*ones(1,dimension));
gradientBasis{2} = h2./(hh2*ones(1,dimension));
gradientBasis{3} = h3./(hh3*ones(1,dimension));


indices1 = find(sum(gradientBasis{1}.^2,2)==0|isnan(sum(gradientBasis{1}.^2,2)));
indices2 = find(sum(gradientBasis{2}.^2,2)==0|isnan(sum(gradientBasis{2}.^2,2)));
indices3 = find(sum(gradientBasis{3}.^2,2)==0|isnan(sum(gradientBasis{3}.^2,2)));

min_norm_grad = min([ ...
  sum(gradientBasis{1}(sum(gradientBasis{1}.^2,2) > 0,:).^2,2); ...
  sum(gradientBasis{2}(sum(gradientBasis{2}.^2,2) > 0,:).^2,2); ...
  sum(gradientBasis{3}(sum(gradientBasis{3}.^2,2) > 0,:).^2,2) ...
  ]);

gradientBasis{1}(indices1,:) = repmat([1 1 1]/min_norm_grad, length(indices1), 1);
gradientBasis{2}(indices2,:) = repmat([1 1 1]/min_norm_grad, length(indices2), 1);
gradientBasis{3}(indices3,:) = repmat([1 1 1]/min_norm_grad, length(indices3), 1);


triangleAreas = sqrt(hh1.*vv)/2;
triangleAreas(isnan(triangleAreas)) = 0;


if dimension == 3
    FaceNormals = cross(w,u);
    FaceNormals = FaceNormals./repmat(sqrt(sum(FaceNormals.^2,2)),1,3);
else
    FaceNormals = [];
end

end
