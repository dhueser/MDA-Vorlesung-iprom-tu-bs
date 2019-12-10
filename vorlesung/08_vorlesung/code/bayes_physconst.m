% by Ehret, date: 2019-12-10
clear all; close all;
% Beispiel: Bayes for a physical constants

% Definition der X-Werte
X = linspace(600,1100,1000);

% Beispiel 
% Prior-Verteilung des Beobachters A
% Normalverteilt mit mu_A = 900 und sigma_a = 20
mu_A = 900;
sigma_A = 20;
DeltaX = X(2)-X(1);
p_A = normpdf(X,mu_A, sigma_A);

% Prior-Verteilung des Beobachters B
% Normalverteilt mit mu_B = 800 und sigma_a = 80
mu_B = 800;
sigma_B = 80;
p_B = normpdf(X,mu_B, sigma_B);

% Plotten der Prior-Verteilungen von A und B
figure(1)
plot(X,p_A,'r-','linewidth',2)
hold on;
plot(X,p_B,'b--','linewidth',2)
grid on;
legend('p_A(\theta) ~ N(900,20^2)','p_B(\theta) ~ N(800,80^2)', ...
'Location','NorthWest')
xlabel('\theta','fontsize',14)
ylabel('p','fontsize',14)

set(gcf, 'PaperUnits', 'centimeters');
x_width=15 ;y_width= 8;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
print '-dpng' prior_A_Prior_B.png

% Annahme der Likelihood (Messung) als Normalverteilung mit l ~ N(850,40^2 ) 
mu_likeli = 850;
sigma_likeli = 40;
likelihood = normpdf(X,mu_likeli, sigma_likeli);

% Plotten der Likelihood
figure(2)
plot(X,likelihood,'k-','linewidth',2)
grid on;
legend('l(\theta|X) \sim N(850,40^2)')
xlabel('\theta','fontsize',14)
ylabel('p','fontsize',14)

set(gcf, 'PaperUnits', 'centimeters');
x_width=15 ;y_width= 8;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
print '-dpng' likelihood_1_Beobachtung.png


% Berechnung der Posterior fuer eine Beobachtung fuer A und B
% Posterior ~ Prior  x Likelihood
p_post_A_1 = p_A .* likelihood;
p_post_B_1 = p_B .* likelihood;

% Normierung des posterior mit Normierungskonstante 
% Nicht ganz korrekt, da eigentlich 
% von -\inf bis +\inf integriert bzw. summiert werden muesste

const_post_A1 = 1/sum(p_post_A_1*DeltaX);
p_post_A_1 = const_post_A1 * p_post_A_1 ;

const_post_B1 = 1/sum(p_post_B_1*DeltaX);
p_post_B_1 = const_post_B1 * p_post_B_1 ;

% Plotten der Posterior
figure(3)
plot(X,p_post_A_1,'r-','linewidth',2)
hold on;
plot(X,p_post_B_1,'b--','linewidth',2)
grid on;
legend('p_A(\theta|X)','p_B(\theta|X) ', ...
'Location','NorthWest')
xlabel('\theta','fontsize',14)
ylabel('p','fontsize',14)

set(gcf, 'PaperUnits', 'centimeters');
x_width=15 ;y_width= 8;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
print '-dpng' posterior_A_posterior_B.png

% Jetzt besteht die Likelihood aus 100 Stichproben
% Bei 100 Stichproben kann die Verteilung schmaler werden
% Wir nehmen eine Normalverteilung N(870, \sigma /sqrt(N)^2) an
N = 100;
sigma_likeli_100 = sigma_likeli / sqrt(N);
mu_likeli_100 = 870;
% Likelihood Funktion fuer 100 Stichproben
likelihood_100 = normpdf(X,mu_likeli_100, sigma_likeli_100);

figure(4)
plot(X,likelihood_100,'k-','linewidth',2)
grid on;
legend('l(\theta|X) ~ N(870,4^2)')
xlabel('\theta','fontsize',14)
ylabel('p','fontsize',14)

% Berechnung der Posterior fuer 100 Beobachtungen fuer A und B
p_post_A_100 = p_A .* likelihood_100;
p_post_B_100 = p_B .* likelihood_100;

% Normierung des posterior mit Normierungskonstante 
% Nicht ganz korrekt, da eigentlich 
% von -\inf bis +\inf integriert bzw. summiert werden muesste
const_post_A_100 = 1/sum(p_post_A_100*DeltaX);
p_post_A_100 = const_post_A_100 * p_post_A_100 ;

const_post_B_100 = 1/sum(p_post_B_100*DeltaX);
p_post_B_100 = const_post_B_100 * p_post_B_100 ;

figure(5)
plot(X,p_post_A_100,'r-','linewidth',2)
hold on;
plot(X,p_post_B_100,'b--','linewidth',2)
grid on;
legend('p_A(\theta|X)','p_B(\theta|X) ', ...
'Location','NorthWest')
xlabel('\theta','fontsize',14)
ylabel('p','fontsize',14)

