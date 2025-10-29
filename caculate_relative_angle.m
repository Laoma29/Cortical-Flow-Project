function [x_reltive_angle,z_reltive_angle]=caculate_relative_angle(raw_optical_flow,Faces，Vertices，VertNormals)

dimension=3;
angle_results_normal= caculate_normal_two_angle(VertNormals);
head_x_angle=angle_results_normal.x_theta;
head_z_angle=angle_results_normal.z_theta;
angle_results= caculate_two_angle(raw_optical_flow); 
OF_x_angle=angle_results.x_theta;
OF_z_angle=angle_results.z_theta;
for i = 1:length(OF_x_angle)
x_reltive_angle(i,:)= angdiff(OF_x_angle(i),head_x_angle(i));
z_reltive_angle(i,:)= angdiff(OF_z_angle(i),head_z_angle(i));
end

function angle_results= caculate_normal_two_angle(flowField) 

flowField=flowField';
  normv=sqrt(sum(flowField.^2,1));
  x=flowField(1,:)./normv;
  y=flowField(2,:)./normv;
  z=flowField(3,:)./normv;
  x_theta = asin(x);%
  z_theta =asin(z);%
  for ii=1:size(flowField,2)
    if y(ii) < 0
      x_theta(ii)= -x_theta(ii) + 2*pi;
    end
  end
  for ii=1:size(flowField,2)
    if y(ii) < 0
      z_theta(ii)= -z_theta(ii) + 2*pi;
    end
  end  
  angle_results.x_theta=x_theta;
  angle_results.z_theta=z_theta;

end
function angle_results= caculate_two_angle(flowField) 

flowField=flowField';
  normv=sqrt(sum(flowField.^2,1));
  x=flowField(1,:)./normv;
  y=flowField(2,:)./normv;
  z=flowField(3,:)./normv;
  x_theta = acos(x);
  z_theta =acos(z);
  for ii=1:size(flowField,2)
    if y(ii) < 0
      x_theta(ii)= -x_theta(ii) + 2*pi;
    end
  end
  for ii=1:size(flowField,2)
    if y(ii) < 0
      z_theta(ii)= -z_theta(ii) + 2*pi;
    end
  end  
  angle_results.x_theta=x_theta;
  angle_results.z_theta=z_theta;

end

end