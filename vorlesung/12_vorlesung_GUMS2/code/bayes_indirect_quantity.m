function bayes_indirect_quantity()
% Hueser, 2017-11-25
%  Kalibrierfaktor
  K_faktor_0 = 0.0925;
  sigma_K = 0.0090;
%  Unsicherheit Messgeraet
  sigma_M = 0.3; % Volt
%
%
  JM = 11;
  sig = 12;
  mue = 500;
  data_M = mue + sig*randn( JM,1);
%
%
  xK = [-4*sigma_K:0.001:4*sigma_K] + K_faktor_0;
  nK = length(xK)
  xM = [-4*sig:0.1:4*sig] + mue;
  nM = length(xM)
  y = [20:0.2:80];
  nY = length(y)
%
  std_M = std(data_M)*1.7;
  std_KM = K_faktor_0 * std_M;
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
  figure(110);
  plot( xK, p_K,'linewidth',2);
  xlabel('Kalibrierfaktor / (Einheit/Volt)', 'fontsize', 14);
  ylabel('pdf', 'fontsize', 14);
  set(gca, 'fontsize', 14);
  set(gcf, 'PaperUnits', 'centimeters');
  x_width=15 ;y_width= 8;
  set(gcf, 'PaperPosition', [0 0 x_width y_width]);
  print(110, 'Kalibrierfaktor.png', '-dpng');
%
% Messungen mit der kleinen Streuung
  p_M_matrix = exp(-0.5 * ( XM/sigma_M ).^2 );
%   Summe ueber alle Messungen
  p_M_1 = sum( p_M_matrix );
  p_M_1 = p_M_1 / sum(p_M_1);
% mit der grossen Streuung
  p_M_matrix = exp(-0.5 * ( XM/std_M ).^2 );
%   Summe ueber alle Messungen
  p_M_2 = sum( p_M_matrix );
  p_M_2 = p_M_2 / sum(p_M_2);
  figure(111);
  plot(xM, p_M_1,'k-', ...
    xM, p_M_2,'r--','linewidth',2);
  xlabel('Rohdaten / Volt', 'fontsize', 14);
  ylabel('pdf', 'fontsize', 14);
  legend('Geraet genau','Geraet ungenau','fontsize',14)
  set(gca, 'fontsize', 14);
  set(gcf, 'PaperUnits', 'centimeters');
  x_width=15 ;y_width= 8;
  set(gcf, 'PaperPosition', [0 0 x_width y_width]);
  print(111, 'Rohdaten.png', '-dpng');
%
% PWahrscheinlichkeit als Funktion der beiden
% direkten Groessen xK und xM
  p_KM_matrix = p_K' * p_M_1;
%  in einen langen Spaltenvektor gebracht
  p_KM_1 = p_KM_matrix(:);
%
  p_KM_matrix = p_K' * p_M_2;
%  in einen langen Spaltenvektor gebracht
  p_KM_2 = p_KM_matrix(:);
%
% Wahrscheinlichkeitsdichteverteilung des Modells
% zur Vereinfachung nehme ich jetzt keine 
% Verteilung um die std_m sondern nur den einen Schaetzwert
  p_y = exp(-0.5 * ( delta_y/std_KM ).^2 );
%
%  posterior
  post_matrix = p_y .* (p_KM_1 * ones(1,nY));
  sumpost = sum(post_matrix(:));
  post_matrix = post_matrix/sumpost;
%
% Marginalverteilung durch Summation ueber xM, xK
  posterior1 = sum(post_matrix);
%
%
%  posterior
  post_matrix = p_y .* (p_KM_2 * ones(1,nY));
  sumpost = sum(post_matrix(:));
  post_matrix = post_matrix/sumpost;
%
% Marginalverteilung durch Summation ueber xM, xK
  posterior2 = sum(post_matrix);
%
%
  p_y2 = exp(-0.5 * ( delta_y/(std_KM*0.003) ).^2 );
%  posterior
  post_matrix = p_y2 .* (p_KM_2 * ones(1,nY));
  sumpost = sum(post_matrix(:));
  post_matrix = post_matrix/sumpost;
%
% Marginalverteilung durch Summation ueber xM, xK
  posterior3 = sum(post_matrix);
%
  [y_KMsort, km_sort] = sort(y_KM);
  p_2 = p_KM_2(km_sort);
  diff_y = diff(y_KMsort);
  mid_p_KM_2 = (p_2(2:nKM)+p_2(1:nKM-1))*0.5;
  iuse = find(diff_y>0);
  meandiff = mean(diff_y(iuse))
  diff_y = diff_y/meandiff;
  pdf = mid_p_KM_2(iuse)./diff_y(iuse);
%
  figure(120);
  plot(diff_y);
  figure(121);
  plot(y_KMsort(iuse),pdf);
  sumpdf = 5.2e6; max(pdf)
  figure(100);
  plot( y_KMsort(iuse), max(posterior1)*pdf/sumpdf, 'g.',...
    y, posterior1,'k-',...
    y, posterior2,'r--',...
    y, posterior3,'b-','linewidth',2);
  xlabel('indirekte Groesse / Einheit', 'fontsize', 14);
  ylabel('pdf', 'fontsize', 14);
  legend('Dirac: Weise/Wuebbler','Geraet genau, Messobj schlecht',...
    'beides ungenau','Geraet ungenau, Modell duenner Peak')
% axis([10 100]);
% set(gca, 'fontsize', 14);
  set(gcf, 'PaperUnits', 'centimeters');
  x_width=15 ;y_width= 8;
  set(gcf, 'PaperPosition', [0 0 x_width y_width]);
  print(100, 'Posterior_all.png',  '-dpng');
end