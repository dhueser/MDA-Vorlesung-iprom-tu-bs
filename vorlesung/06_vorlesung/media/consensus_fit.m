function consensus_fit()
  x_A = [2.0,1.0,1.5,1.8,1.2,1.7];
  x_B = [16.3,16.8];

% jeweiliger Stichprobenumfang
  N_A = length(x_A);
  N_B = length(x_B);

% Mittelwerte
  bar_x_A = mean(x_A);
  bar_x_B = mean(x_B);
  fprintf(stdout,'Mittelwerte: %1.2f und %1.2f\n', ...
    bar_x_A, bar_x_B);

% Varianzen der Mittelwerte
  Var_A_mittel = var(x_A)/N_A;
  Var_B_mittel = var(x_B)/N_B;
  fprintf(stdout,'Varianzen der Mittelwerte innherhalb jeder Stichprobe:\n %1.3f und %1.3f\n', ...
    Var_A_mittel, Var_B_mittel);

% Berechnung des Birge ratio fuer den Startwert v = 0.0
  v = 0.0;
  [RB, x_ref, dv] = eval_Birge_dv([Var_A_mittel, Var_B_mittel], ...
       [bar_x_A, bar_x_B], v);
% 
  max_num_iterations = 50;
  k = 0;
  while (RB > 1) && (k < max_num_iterations)
    fprintf(stdout,'x_ref: %1.3f, RB: %1.2f - dv: %e ==> v = %1.3f\n', ...
       x_ref, RB, dv, v);
    v = v + dv;
    [RB, x_ref, dv] = eval_Birge_dv([Var_A_mittel, Var_B_mittel], ...
       [bar_x_A, bar_x_B], v);
    k = k + 1;
  end
  fprintf(stdout,'x_ref: %1.3f, RB: %1.2f - dv: %e ==> v = %1.3f\n', ...
       x_ref, RB, dv, v);
  fprintf(stdout,'Standardabweichung zwischen den Partnern:\n');
  fprintf(stdout,'s_b = %1.1f nach %d Iterationen\n', ...
    sqrt(v), k);
%
  var_pool = eval_poolvar([Var_A_mittel, Var_B_mittel], ...
       [bar_x_A, bar_x_B], v);
  fprintf(stdout,'var_pool = %1.3f == std_pool = %1.1f\n', ...
    var_pool, sqrt(var_pool));
end
%
%
function [RB, x_ref, dv] = eval_Birge_dv(var_barx, bar_x, var_b)
  n = length(bar_x);
  weight = 1./(var_barx + var_b);
  x_ref = (1/(n-1)) * sum(weight.*bar_x) / sum(weight);
  RBsq = (1/(n-1)) * sum(weight .* (bar_x-x_ref).^2);
  RB = sqrt(RBsq);
%
  Q = RBsq - 1;
  dQ_dv = (-1/(n-1)) * sum(weight.^2 .* (bar_x-x_ref).^2);
  dv = -Q / dQ_dv;
end
%
%
function var_pool = eval_poolvar(var_barx, bar_x, var_b)
  weight = 1./(var_barx + var_b);
  var_pool = 1 / sum(weight);
end