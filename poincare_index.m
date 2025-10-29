function index = poincare_index(flowField)
  theta = myangle(flowField);
  difftheta = diffangle([theta(2),theta(3),theta(1)], theta);
  index = sum(difftheta)/(2*pi);
end
