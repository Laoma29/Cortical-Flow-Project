function angle_results= caculate_two_angle(flowField) 

flowField=flowField';
  normv=sqrt(sum(flowField.^2,1));
  x=flowField(1,:)./normv;%x 方向
  y=flowField(2,:)./normv;%y 方向
  z=flowField(3,:)./normv;%y 方向
  x_theta = acos(x);%x方向角度
  z_theta = acos(z);%z方向角度 acos
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