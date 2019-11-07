close all; clear all;
% Eingabe der Messwertepaare
input = [...
1, 0.4;...
2, 0.55;...
3, 0.7;...
4, 0.75; ...
5, 0.8];
m = length(input); % Anzahl der Messpunkte
x_Werte = input(:,1);
X = [ones(m,1), x_Werte, x_Werte.^2,x_Werte.^3];
Y = input(:,2);
% Aufstellen der X-Matrix % Aufstellen der Y-Matrix
% degree of freedom
df = 5-4;
% Berechung der Schätzwerte \thetaHat;
% y-Abschnitt = \thetaHat(1); Steigung = \thetaHat(2)
thetaHat = (X' * X) \ (X' * Y);
% Anzeige der Schätzwerte
disp(thetaHat);
% Plotten der Fitfunktion
xx = linspace(0,6,100);
yy = thetaHat(1) + thetaHat(2)*xx+ thetaHat(3)*xx.^2+thetaHat(4)*xx.^3;
figure(10)
plot(xx, yy,'linewidth',1.5)
% Plotten der Datenwerte
hold on
plot(input(:,1), input(:,2), 'or','linewidth',2, 'MarkerSize',6)
hold off
grid on;
xlabel('Werte X_j','fontsize',14)
ylabel('Werte Y_j','fontsize',14)
set(gcf, 'PaperUnits', 'centimeters');
x_width=15 ;y_width= 9;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
print '-dpng' kubischerFit.png
% Qualitätsmaß Q
Q = (Y-X*thetaHat)' * (Y-X*thetaHat)
% Varianzen der Residuen
Var_epsilon = 1/df *(Y-X*thetaHat)' * (Y-X*thetaHat)
% Wuzel aus Var_epsilon
RMSE = sqrt(Var_epsilon)
% Kovarianzen
Cov = Var_epsilon * inv(X'*X)
