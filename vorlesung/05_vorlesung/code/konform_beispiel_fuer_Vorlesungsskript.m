function konform_beispiel()
  step = 0.001;
  y0 = 1500;
  s0 = 0.12;
  T_L = 1499.80;
  T_U = 1500.20;
  A_L = 1499.82;
  A_U = 1500.18;
  mininf = T_L - 3*s0;
  maxinf = T_U + 3*s0;

  y = [T_L:step:T_U];
  prior = step * exp(-0.5*((y - y0) / s0).^2) / (s0 * sqrt(2*pi));
  % conformance probability p_c
  p_c = sum(prior)

% inspection device - instrument uncertainty
  s_insp = 0.04;
  yall = [mininf:step:maxinf];
  resistor_distri = exp(-0.5*((yall - y0) / s0).^2) / (s0 * sqrt(2*pi));
  p_TL = exp(-0.5*((T_L - y0) / s0).^2) / (s0 * sqrt(2*pi));
  p_TU = exp(-0.5*((T_U - y0) / s0).^2) / (s0 * sqrt(2*pi));
  p_AL = exp(-0.5*((A_L - y0) / s0).^2) / (s0 * sqrt(2*pi));
  p_AU = exp(-0.5*((A_U - y0) / s0).^2) / (s0 * sqrt(2*pi));

  x1 = T_U + 0.03;
  x = [x1-5*s_insp:step:x1+5*s_insp];
  p_insp = (exp(-0.5*((x - x1) / s_insp).^2) / (s_insp * sqrt(2*pi))) * ...
     (exp(-0.5*((x1 - y0) / s0).^2) / (s0 * sqrt(2*pi)));
  hlp = (1 / (s_insp * sqrt(2*pi))) * ...
     (exp(-0.5*((x1 - y0) / s0).^2) / (s0 * sqrt(2*pi)));
  hlpAU = (exp(-0.5*((A_U - x1) / s_insp).^2) / (s_insp * sqrt(2*pi))) * ...
     (exp(-0.5*((x1 - y0) / s0).^2) / (s0 * sqrt(2*pi)));

% producer risk R_P:
  x_L = [mininf:step:A_L]'; 
  x_U = [A_U:step:maxinf]';
  n_y = length(y);
  n_xL = length(x_L);
  n_xU = length(x_U);
  likeli_L = exp(-0.5 * ((x_L * ones(1,n_y) - ones(n_xL,1) * y) / s_insp).^2) ...
             * step / (s_insp * sqrt(2*pi));
  likeli_U = exp(-0.5 * ((x_U * ones(1,n_y) - ones(n_xU,1) * y) / s_insp).^2) ...
             * step / (s_insp * sqrt(2*pi));

  posterior_L = (ones(n_xL,1)*prior) .* likeli_L;
  posterior_U = (ones(n_xU,1)*prior) .* likeli_U;
  sum_L = sum(posterior_L(:));
  sum_U = sum(posterior_U(:));
  R_P = sum_L + sum_U

% consumer risk R_C
  x = [A_L:step:A_U]'; 
  y_L = [mininf:step:T_L];
  y_U = [T_U:step:maxinf];
  n_x = length(x);
  n_yL = length(y_L);
  n_yU = length(y_U);
  prior_L = step * ones(n_x,1) * exp(-0.5*((y_L - y0) / s0).^2) / (s0 * sqrt(2*pi));
  prior_U = step * ones(n_x,1) * exp(-0.5*((y_U - y0) / s0).^2) / (s0 * sqrt(2*pi));
  likeli_L = exp(-0.5 * ((x * ones(1,n_yL) - ones(n_x,1) * y_L) / s_insp).^2) ...
             * step / (s_insp * sqrt(2*pi));
  likeli_U = exp(-0.5 * ((x * ones(1,n_yU) - ones(n_x,1) * y_U) / s_insp).^2) ...
             * step / (s_insp * sqrt(2*pi));
  posterior_L = prior_L .* likeli_L;
  posterior_U = prior_U .* likeli_U;
  sum_L = sum(posterior_L(:));
  sum_U = sum(posterior_U(:));
  R_C = sum_L + sum_U
end
