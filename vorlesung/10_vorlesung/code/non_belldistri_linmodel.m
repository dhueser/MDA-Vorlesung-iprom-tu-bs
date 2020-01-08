function non_belldistri_linmodel (sigma=0.05)
% Hueser, 2019-12-19 
  
% sigma = 0.05, d.h. der Modellparanmter y wird hier mit gaussfoermigem 
% Modellprior, aber sehr schmal, fast dirac'sche delta-Distribution angenommen.
  
  close all
  
% Falls das Statistikpaket noch nicht geladen ist, kann es wie folgt geladen werden:
% Entsprechend die folgende Zeile ein- oder auskommentieren
  pkg load statistics

  boxflag = 0; % 0: U-Verteilung; 1: Rechteck-Verteilung;
  
  pltflag1 = 1;
  pltflag1a = 0;
  pltflag2 = 1;
  pltflag3 = 1;
   
  a = 1.0;
  b = 1.0;
  x1_quer = 3.0;
% Spanne
  s1 = 2.0;
%
  x2_quer = 3.5;
% Standardabweichung
  s2 = 0.5;
%
% Modellgleichung y = a x1 + b x2
%
% Zeilenvektor
  x1 = [-s1:0.05:s1] + x1_quer;
  if boxflag
    p1 = ones(size(x1));
  else
    p1 = asin((x1 - x1_quer)/s1).^2;
  end
  p1 = p1/sum(p1);
  figure(1);
  plot(x1, p1, 'k-');
  grid on;
  xlabel('direkte Messgroesse x1', 'fontsize', 14);
  ylabel('Wahrscheinlichkeitsdichte p(x1)', 'fontsize', 14);
  set(gca, 'fontsize', 12);
  if pltflag1
    print(1,'direkte_x1_box.png', '-dpng');
  end
%
% Spaltenvektor
  x2 = [-3*s2:0.05:3*s2]' + x2_quer;
  p2 = normpdf(x2, x2_quer, s2);
  p2chk = exp(-0.5 * ((x2 - x2_quer)/s2).^2) / (sqrt(2*pi) * s2);
  figure(2);
  plot(x2, p2, 'k-', x2, p2chk, 'r--');
  grid on;
  xlabel('direkte Messgroesse x2', 'fontsize', 14);
  ylabel('Wahrscheinlichkeitsdichte p(x2)', 'fontsize', 14);
  set(gca, 'fontsize', 12);
  if pltflag2
    print(2,'direkte_x2.png', '-dpng');
  end
%
% Spalte mal Zeile liefert Matrix
  pgesamt = p2*p1;
  pgesamt = pgesamt / sum(pgesamt(:));
%
  figure(3);
  imagesc(x1, x2, pgesamt);
  colorbar();
  xlabel('direkte Messgroesse x1', 'fontsize', 14);
  ylabel('direkte Messgroesse x2', 'fontsize', 14);
  title('Wahrscheinlichkeitsdichteverteilung p(x1, x2)', 'fontsize', 14);
  set(gca, 'fontsize', 12);
  if pltflag1
    print(3,'direkte_x1_x2_box.png', '-dpng');
  end
%
% ----
%  x1 und x2 in Matrizen fuellen
  n2_r = length(x2)
  n1_c = length(x1)
  X_1 = ones(n2_r,1) * x1;
  X_2 = x2 * ones(1, n1_c);
  Y_matrix = a*X_1 + b*X_2;
%
% nun als 1D arrays
  y = Y_matrix(:);
  p_y = pgesamt(:);
  figure(4);
  plot(y, p_y);
  axis([2 11]);
  title('scharfer Modellprior - GUM S1', 'fontsize', 14);
  xlabel('indirekte Messgroesse y', 'fontsize', 14);
  ylabel('Wahrscheinlichkeitsdichte', 'fontsize', 14);
  set(gca, 'fontsize', 12);
  if pltflag1a
    print(4,'indirekte_y_GUMS1.png', '-dpng');
  end
%
% sortieren, so dass y monoton steigend ist
  [ys, isort] = sort(y);
  p_sort = p_y(isort);
  figure(5);
  plot(ys, p_sort, 'x-');
  axis([2 11]);
  grid on;
  title('scharfer Modellprior - GUM S1', 'fontsize', 14);
  xlabel('indirekte Messgroesse y', 'fontsize', 14);
  ylabel('Wahrscheinlichkeitsdichte', 'fontsize', 14);
  set(gca, 'fontsize', 12);
  if pltflag1
    print(5,'indirekte_ysorted_GUMS1_box.png', '-dpng');
  end

%
% mit gaussfoermigem Modellprior, aber sehr schmal, fast dirac'sche delta-Distribution
%
  dy = 0.01;
  y_model = [2:dy:11];
  n_model = length(y_model);
%  sigma = 0.5;
  p_mod = zeros(1, n_model);
  for k = 1:n_model
%
% Summationen entsprechen der Integration ueber x1 und x2
% zum Berechnen der Marginalverteilung
    p_mod(k) = sum(sum(exp(-0.5 * ((y_model(k) - Y_matrix)/sigma).^2) .* pgesamt));
  end
  Pkum = cumsum(p_mod);
  psum = Pkum(n_model);
  Pkum = Pkum/psum;
  p_mod = p_mod/psum;

%
% Signifikanzniveau
  alpha = 0.05;
%
% Quantile
  [left, ileft] = min(abs(Pkum - alpha/2))
  [right, iright] = min(abs(Pkum - (1-alpha/2)))
%
  erw_y = sum(y_model .* p_mod);
  std_y = sqrt(sum((y_model - erw_y).^2 .* p_mod));
%
% Vergleich mit Fortpflanzungsgesetz, das eigentlich nur fuer
% glockenfoermige Verteilungen gilt, bei denen die Wurzel des zweiten statistischen
% Moments, also die Varianz die Breite der Verteilung definiert, bzgl der sie
% sich normieren laesst.
  y_quer = a*x1_quer + b*x2_quer;
  if boxflag
    u_y = sqrt((a * s1)^2/3 + (b * s2)^2);
  else
    u_y = sqrt((a * s1)^2/2 + (b * s2)^2);
  end
%
  figure(6);
  sigma_str = num2str(sigma,'%1.3f');
  plot(y_model, p_mod, 'kx-', ...
      [y_quer-1.96*u_y, y_quer-1.96*u_y], [0, p_mod(iright)], 'b-.', 'linewidth', 2.5, ...
      [y_quer+1.96*u_y, y_quer+1.96*u_y], [0, p_mod(iright)], 'b-.', 'linewidth', 2.5, ...
      [y_model(ileft), y_model(ileft)], [0, p_mod(iright)], 'r--', 'linewidth', 3, ...
      [y_model(iright), y_model(iright)], [0, p_mod(iright)], 'r--', 'linewidth', 3);
  grid on;
  axis([2 11]);
  title(['normalverteilter Modellprior mit \sigma = ', sigma_str], 'fontsize', 14);
  xlabel('indirekte Messgroesse y', 'fontsize', 14);
  ylabel('Wahrscheinlichkeitsdichte', 'fontsize', 14);
  set(gca, 'fontsize', 12);
  if pltflag3
    print(6,['indirekte_y_modellGauss_sigma', strrep(sigma_str,'.','p'),'_box.png'], '-dpng');
  end

%
  figure(7);
  plot(y_model, Pkum, 'ko-', ...
    [y_model(ileft), y_model(ileft)], [0, Pkum(ileft)], 'r--', 'linewidth', 3, ...
    [y_model(iright), y_model(iright)], [0, Pkum(iright)], 'r--', 'linewidth', 3, ...
    [y_model(1), y_model(ileft)], [Pkum(ileft), Pkum(ileft)], 'r--', 'linewidth', 3, ...
    [y_model(1), y_model(iright)], [Pkum(iright), Pkum(iright)], 'r--', 'linewidth', 3);
  grid on;
  axis([2 11]);
  title(['normalverteilter Modellprior mit \sigma = ', sigma_str], 'fontsize', 14);
  xlabel('indirekte Messgroesse y', 'fontsize', 14);
  ylabel('kumulierte Wahrscheinlichkeit', 'fontsize', 14);
  set(gca, 'fontsize', 12);
  if pltflag3
    print(7,['Pkum_indirekte_y_Gauss_sigma', strrep(sigma_str,'.','p'),'_box.png'], '-dpng');
  end
%
  printf('Modellprior mit sigma = %1.3f\n', sigma);
  printf('Erwartungswerte:\n');
  printf('aus Erw von x1 und x2:  %1.2f\n', y_quer);
  printf('aus Verteilung p_model: %1.2f\n', erw_y);
  printf('Standardabweichungen als Wurzel der Varianz:\n');
  printf('aus Fortpflanzung klassisch: %1.2f\n', u_y);
  printf('aus Verteilung p_model:      %1.2f\n', std_y);
  printf('Ueberdeckungsintervall fuer 95 Prozent\n');
  printf('aus Fortpflanzung klassisch mit Quantil aus Gaussvert 1.96:  [%1.2f, %1.2f]\n',
    y_quer - 1.96*u_y, y_quer + 1.96*u_y);
  printf('aus Verteilung p_model, ueber die Umkehrfkt der kumulativen: [%1.2f, %1.2f]\n',
    y_model(ileft), y_model(iright));
end