function new_vector_global=Projection_new_vector(flowField_standard,tangentPlaneBasis)
all_x_vector=tangentPlaneBasis(:,:,1);
all_y_vector=tangentPlaneBasis(:,:,2);
all_z_vector=tangentPlaneBasis(:,:,3);
j=1;
raw_OF = mean(flowField_standard,3);
% for j = 1:size(raw_OF,3)
                syms x y z
parfor i = 1:size(raw_OF,1)
            raw_OF_1 = squeeze(raw_OF(i,1));
    raw_OF_2 = squeeze(raw_OF(i,2));
    raw_OF_3 = squeeze(raw_OF(i,3));
[x1,y1,z1]=solve(x*all_x_vector(i,1)+y*all_y_vector(i,1)+z*all_z_vector(i,1)==raw_OF_1,...
    x*all_x_vector(i,2)+y*all_y_vector(i,2)+z*all_z_vector(i,2)==raw_OF_2,...
    x*all_x_vector(i,3)+y*all_y_vector(i,3)+z*all_z_vector(i,3)==raw_OF_3);
x_value =subs(x1);
x_value=double(x_value);
y_value =subs(y1);
y_value=double(y_value);
z_value =subs(z1);
z_value=double(z_value);
new_vector(i,:)=[x_value y_value z_value];
end