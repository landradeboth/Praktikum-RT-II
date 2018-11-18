clear all;
close all;
clc;

% syms k T;
% A = [0 1; -k/T -2*sqrt(k/T)];
% 
% [V,D] = eig(A);

ks = -6.4;
Tm = 0.27;
ki = -7.27;
freq_rad = sqrt(10/Tm);

%% Zweipunktglied
A_Zwei = (Tm*ks*ki)/(5*pi*(Tm+0.1));

%% Dreipunktglied
C = ((Tm*ks*ki)/(5*pi*(Tm + 0.1)))^2;
A1 = sqrt((C + sqrt(C^2 - C))/2);
A2 = sqrt((C - sqrt(C^2 - C))/2);


H = tf([0.5*ks*ki],[0.1*Tm (Tm + 0.1) 1 0]);
w = linspace(0,150,1000);
nyquist(H,w);
axis([-9 1 -20 20]);