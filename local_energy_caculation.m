function [Local_dEnergy,Global_kinetic_energy]=local_energy_caculation(Vertices，Faces,VertNormals，flowField)
tic;

dimension=3;
[syed, triangleAreas] = geometry_tesselation(Faces, Vertices, dimension);
dEnergy = zeros(1, size(flowField,3));
        v12 = squeeze(sum((flowField(Faces(:,1),:,:)+flowField(Faces(:,2),:,:)).^2,2) / 4);
        v23 =squeeze( sum((flowField(Faces(:,2),:,:)+flowField(Faces(:,3),:,:)).^2,2) / 4);
        v13 = squeeze(sum((flowField(Faces(:,1),:,:)+flowField(Faces(:,3),:,:)).^2,2) / 4);
        Local_dEnergy1 = (triangleAreas.*(v12+v23+v13));%%计算所有三角形的动能
        Global_kinetic_energy=sum(Local_dEnergy1);
    for i =1:size(Vertices,1)
   Local_dEnergy(i,:)=  (sum(Local_dEnergy1(Faces(:,1)==i,:),1)+sum(Local_dEnergy1(Faces(:,2)==i,:),1)+sum(Local_dEnergy1(Faces(:,3)==i,:),1))/sum(sum(Faces==i));%

    end
    
    Local_dEnergy = sqrt(Local_dEnergy); % Get square root for easier visualization
t=toc;
disp(t);
end