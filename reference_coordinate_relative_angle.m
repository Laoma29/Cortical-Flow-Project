function relative_angle=reference_coordinate_relative_angle(raw_OF,FV)

Faces=FV.Faces;
Vertices= FV.Vertices;
VertNormals=FV.VertNormals;
dimension=3;
[gradientBasis, triangleAreas, FaceNormals] = ...
  geometry_tesselation(Faces, Vertices, dimension);
%%
[regularizerOld, gradientBasis, tangentPlaneBasisCell, tangentPlaneBasis, ...
  triangleAreas, FaceNormals, sparseIndex1, sparseIndex2]= ...
  regularizing_matrix(Faces, Vertices, VertNormals, dimension); 
x_vector=tangentPlaneBasis(:,:,1);
y_vector=tangentPlaneBasis(:,:,2);
z_vector=-VertNormals;
for i = 1:size(raw_OF,1)
    syms x y z
[x,y,z]=solve(x*x_vector(i,1)+y*y_vector(i,1)+z*z_vector(i,1)==raw_OF(i,1),...
    x*x_vector(i,2)+y*y_vector(i,2)+z*z_vector(i,2)==raw_OF(i,2),...
    x*x_vector(i,3)+y*y_vector(i,3)+z*z_vector(i,3)==raw_OF(i,3));
x_value =subs(x);
x_value=double(x_value);
y_value =subs(y);
y_value=double(y_value);
z_value =subs(z);
z_value=double(z_value);
new_vector(i,:)=[x_value y_value z_value];
end
for i = 1:length(new_vector)
relative_angle(i,:)=acos(dot(new_vector(i,:),x_vector(i,:))/(norm(new_vector(i,:))*norm(x_vector(i,:))));
if new_vector(i,2)<0
relative_angle(i,:)=2*pi-relative_angle(i,:);
end
end


end