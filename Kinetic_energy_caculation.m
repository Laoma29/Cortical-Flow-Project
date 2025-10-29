function [Global_dEnergy, Local_dEnergy,dspeed]=Kinetic_energy_caculation(Faces,Vertices，VertNormals，flowField)

dimension=3;
[syed, triangleAreas] = geometry_tesselation(Faces, Vertices, dimension);
dEnergy = zeros(1, size(flowField,3));
dspeed  = zeros(size(flowField,1), size(flowField,3));
% dEnergy2 = zeros(1, size(flowField,3));
for m = 1:size(flowField,3)
    Global_dEnergy(:,m)=sum(sum(triangleAreas.*sum((flowField(Faces(:,1),:,m)+flowField(Faces(:,2),:,m)+flowField(Faces(:,3),:,m)).^2)));
    Local_dEnergy1(:,m)=sum(triangleAreas.*sum((flowField(Faces(:,1),:,m)+flowField(Faces(:,2),:,m)+flowField(Faces(:,3),:,m)).^2).^2,2);

    for i =1:size(Vertices)
        Local_dEnergy(i,m)=  (sum(sum(Local_dEnergy1(Faces(:,1)==i,m)))+sum(sum(Local_dEnergy1(Faces(:,2)==i,m)))+sum(sum(Local_dEnergy1(Faces(:,3))==i,m)))/sum(sum(Faces==i));
    end

    dspeed(:,m)=sqrt(flowField(:,1,m).^2+flowField(:,2,m).^2+flowField(:,3,m).^2);

    
end
Global_dEnergy = sqrt(Global_dEnergy); % Get square root for easier visualization
Local_dEnergy = sqrt(Local_dEnergy); % Get square root for easier visualization
dspeed = dspeed;