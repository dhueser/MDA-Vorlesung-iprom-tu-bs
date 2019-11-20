function example_t()
  J1 = 50;
  mu1 = 42;
  sigma1 = 5;
  samples = [];
  for kappa = 1:5
    samples = [samples  (mu1 + sigma1*randn(J1,1))];
  end
  figure(3000); hold on;
  [J1, nS] = size(samples);
  for kappa = 1:5
    [P, H, xsort, xbar, s] = cumulatives(samples(:,kappa));
    plot( xsort, H, 'k.-');
    plot( xsort, P, 'r-');
  end
  xlabel('X', 'fontsize', 14);
  ylabel('Wahrscheinlichkeit', 'fontsize', 12);
  title(['Umfang jeweilige Stichproben: ' num2str(J1)], 'fontsize', 14);
  grid on;
  set(gca, 'fontsize', 14, 'linewidth', 2);
  hold off;
end
function [P, H, xsort, xbar, s] = cumulatives(x)
  xbar = mean(x);
  s = std(x);
  J1 = length(x);
  [xsort, isort] = sort(x, 'ascend');
  H = [1:J1]'/J1;
  P = normcdf( xsort, xbar, s);
end