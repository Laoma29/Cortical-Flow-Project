function theta = caculate_angle(flowField) 
flowField=flowField';
  normv=sqrt(sum(flowField.^2,1));
  c=flowField(1,:)./normv;
  s=flowField(2,:)./normv;
  theta = acos(c);
  for ii=1:size(flowField,2)
    if s(ii) < 0
      theta(ii)= -theta(ii) + 2*pi;
    end
  end
end