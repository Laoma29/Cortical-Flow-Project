function angle_results= single_vector_angle_caculation(vector) 
vector=vector';

  normv=sqrt(sum(vector.^2,1));
  x=vector(1,:)./normv;%
  y=vector(2,:)./normv;%
  z=vector(3,:)./normv;%
  x_theta = acos(x);%
  z_theta =acos(z);%
  for ii=1:size(vector,2)
    if y(ii) < 0
      x_theta(ii)= -x_theta(ii) + 2*pi;
    end
  end
  for ii=1:size(vector,2)
    if y(ii) < 0
      z_theta(ii)= -z_theta(ii) + 2*pi;
    end
  end  
  angle_results.x_theta=x_theta;
  angle_results.z_theta=z_theta;

end