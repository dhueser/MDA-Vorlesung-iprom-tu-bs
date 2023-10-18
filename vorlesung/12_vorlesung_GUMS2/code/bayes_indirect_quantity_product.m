function bayes_indirect_quantity_product()
%  Kalibrierfaktor
  K_faktor_0 = 0.0925;
  sigma_K = 0.0090;
%
%
  JM = 9;
  sig = 12;
  mue = 500;
  data_M = [479.58; 526.47; 516.77; 522.01; 506.61; 497.99; 481.71; 484.90; 491.41];
%
%
  xK = [-4*sigma_K:0.005:4*sigma_K] + K_faktor_0;
  nK = length(xK)
  xM = [-4*sig:0.05:4*sig] + mue;
  nM = length(xM)
  y = [25:0.005:65];
  nY = length(y)
%
  std_M = std(data_M);
  std_KM = K_faktor_0 * std_M;
  printf('s_M = %1.4f, s_KM = %1.4f\n', std_M, std_KM);
%
% Messungen in Volt
  XM = data_M * ones(1,nM) - ones(JM,1) * xM;
%
% Das Modell: y = f(xK, xM)
% fuer jedes xK und jedes xM kombiniert
% y_{KM,i,j} = xK_i * xM_j fuer alle i=1,..,nK und j=1,..,nM
  y_KM_matrix = xK' * xM;
%  in einen langen Spaltenvektor gebracht
  y_KM = y_KM_matrix(:);
  nKM = length(y_KM);
  delta_y = y_KM * ones(1,nY) - ones(nKM, 1)*y;
%
% direkte Groessen
% Kalibrierfaktor
  p_K = exp(-0.5 * ( (xK - K_faktor_0)/sigma_K ).^2 );
  p_K = p_K / sum(p_K);
% Messungen: Likelihood
  p_M_matrix = exp(-0.5 * ( XM/std_M ).^2 );
% Summe ueber alle Messungen
  p_M_2 = sum( p_M_matrix );
  p_M_2 = p_M_2 / sum(p_M_2);
%
  p_KM_matrix = p_K' * p_M_2;
%  in einen langen Spaltenvektor gebracht
  p_KM_2 = p_KM_matrix(:);
%
% Modellprior
  p_y2 = exp(-0.5 * ( delta_y/(std_KM*0.003) ).^2 );
% Posterior
  post_matrix = p_y2 .* (p_KM_2 * ones(1,nY));
% Normierung des Posteriors
  sumpost = sum(post_matrix(:));
  post_matrix = post_matrix/sumpost;
%
% Marginalverteilung durch Summation ueber xM, xK
  posterior3 = sum(post_matrix);
%
% kumulative Marginale-PDF also die CDF
  cdf3 = cumsum(posterior3);
%
  figure(200);
  plot( y, cdf3, 'b-', 'linewidth', 2);
  xlabel('indirekte Groesse / Einheit', 'fontsize', 14);
  ylabel('cdf', 'fontsize', 14);
  set(gca, 'fontsize', 12);
%
% Ergebnis
  y_bar = sum( y.*posterior3 );
  [dmy, imin] = min( abs(cdf3-0.025) )
  y1 = y(imin);
  [dmy, imin] = min( abs(cdf3-0.975) )
  y2 = y(imin);
  printf('Credible interval = [%1.4f, %1.4f]\n', y1, y2);
  printf('y_bar = %1.4f + (%1.4f) + (%1.4f)\n', y_bar, y1-y_bar, y2-y_bar);
% ---
% Vergleich mit Fortpflanzungsgesetz
% ---
  xMbar = mean(data_M);
  printf('mean(X_M): %1.4f\n', xMbar);
  printf('s_M = %1.4f, s_KM = %1.4f\n', std_M, std_KM);
  u = sqrt((K_faktor_0*std_M)^2 + (xMbar*sigma_K)^2);
  printf('u(Y) = %1.4f\n',u);
  hlp = (K_faktor_0*std_M)^4/8 + (xMbar*sigma_K)^4/45;
  nu_eff = u^4/hlp
  t = tinv(0.975,floor(nu_eff))
  printf('U = %1.5f \n', t*u);
  printf('ymean = %1.3f, [%1.3f,%1.3f]\n', ...
    K_faktor_0*xMbar, K_faktor_0*xMbar-t*u,K_faktor_0*xMbar+t*u);
end