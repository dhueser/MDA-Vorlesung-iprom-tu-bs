close all; clear all;
% Eingabe der Messwertepaare
input = [...                 
1, 0.4;...                
2, 0.55;...
3, 0.7;...
4, 0.75; ...
5, 0.8];
m = length(input);             % Anzahl der Messpunkte
X = [ones(m,1), input(:,1)];   % Aufstellen der X-Matrix
Y = input(:,2);                % Aufstellen der Y-Matrix
% degree of freedom 
df = 5-2;
% Berechung der Schaetzwerte \thetaHat; 
% y-Abschnitt = \thetaHat(1); Steigung = \thetaHat(2)
thetaHat = (X' * X) \ (X' * Y);    
% Anzeige der Schaetzwerte
disp(thetaHat);
% Plotten der Regressionsgeraden
xx = linspace(0,6,100);
yy = thetaHat(1) + thetaHat(2)*xx;
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
print '-dpng' LinRegression.png
