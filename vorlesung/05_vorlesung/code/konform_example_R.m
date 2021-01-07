function konform_example_R()
  step = 0.0002;
  perc = '%';

  gaussian = @(x, mu, sigma) exp(-0.5*((x - mu) / sigma).^2) / (sigma * sqrt(2*pi));

% ----
  y0 = 1500;
  s0 = 0.12;
  T_L = 1499.80;
  T_U = 1500.20;
  A_L = 1499.82;
  A_U = 1500.18;
  mininf = T_L - 3*s0;
  maxinf = T_U + 3*s0;

  y = [A_L:step:A_U];
  prior = step * gaussian(y, y0, s0);
  p_c = sum(prior);
  printf('spezifisches Konsumentenrisiko bzgl A ..................: %3.0f %c\n', 100 - 100*p_c, perc);


  y = [T_L:step:T_U];
  prior = step * gaussian(y, y0, s0);
  p_c = sum(prior);
  printf('spezifisches Konsumentenrisiko bzgl T Gl (21) in JCGM106: %3.0f %c\n', 100 - 100*p_c, perc);

% inspection device - instrument uncertainty
  s_insp = 0.04;

% ------
% global producer risk:
  x_L = [mininf:step:A_L]'; %'
  x_U = [A_U:step:maxinf]'; %'

  n_y = length(y);
  n_xL = length(x_L);
  n_xU = length(x_U);

  px_giveny_L = step * gaussian(x_L * ones(1,n_y), ones(n_xL,1) * y, s_insp);
  px_giveny_U = step * gaussian(x_U * ones(1,n_y), ones(n_xU,1) * y, s_insp);

  posterior_L = (ones(n_xL,1)*prior) .* px_giveny_L;
  posterior_U = (ones(n_xU,1)*prior) .* px_giveny_U;
  R_p = sum(posterior_L(:)) + sum(posterior_U(:));
  printf('globales Produzentenrisiko: %3.0f %c\n', 100*R_p, perc);

%-------
% global consumer risk

  x = [A_L:step:A_U]'; %'
  y_L = [mininf:step:T_L];
  y_U = [T_U:step:maxinf];

  n_x = length(x);
  n_yL = length(y_L);
  n_yU = length(y_U);
  prior_L = step * ones(n_x,1) * gaussian(y_L, y0, s0);
  prior_U = step * ones(n_x,1) * gaussian(y_U, y0, s0);
  px_giveny_L = step * gaussian(x * ones(1,n_yL), ones(n_x,1) * y_L, s_insp);
  px_giveny_U = step * gaussian(x * ones(1,n_yU), ones(n_x,1) * y_U, s_insp);
  posterior_L = prior_L .* px_giveny_L;
  posterior_U = prior_U .* px_giveny_U;
  R_c = sum(posterior_L(:)) + sum(posterior_U(:));
  printf('globales Konsumentenrisiko: %3.0f %c\n', 100*R_c, perc);
end
