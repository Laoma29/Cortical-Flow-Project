function theta = diffangle(theta2, theta1)
  theta = theta2 - theta1;
  for ii=1:length(theta)
    if theta(ii) < -pi
      theta(ii) = theta(ii) + 2*pi;
    elseif theta(ii) > pi
        theta(ii) = theta(ii) - 2*pi;
    end
  end
end