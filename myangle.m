
% gives the good angle of a vector (between 0 and 2pi)
function theta = myangle(flowField) 
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
