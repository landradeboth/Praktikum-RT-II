clear all;
close all;
clc;

syms m_K m_L l J R K_E K_M K_U K_red r x_max v_max U_PC_max g;

syms m_K_schlange;

%% Aufgabe 9

A = [0                                         1                                        0 0; ...
     0   -(K_M*K_E/(m_K_schlange*R))*(K_red/r)^2                    -(m_L/m_K_schlange)*g 0; ...
     0                                         0                                        0 1; ...
     0 -(K_M*K_E/(m_K_schlange*R*l))*(K_red/r)^2 -((m_L+m_K_schlange)/(m_K_schlange*l))*g 0];

b = [0; (K_M*K_red*K_U)/(R*r*m_K_schlange); 0; (K_M*K_red*K_U)/(l*m_K_schlange*R*r)];

C = [1 0 0 0; 0 0 1 0];

%% Aufgabe 10
% Aufgabe 10a
C_a = [1 0 -l 0];

Beo_Matrix1 = [C_a; C_a*A; C_a*A*A; C_a*A*A*A];

det_Beo_Matrix1 = det(Beo_Matrix1);

rang_Beo_Matrix1 = rank(Beo_Matrix1);

C_a = double([1 0 -1 0]);

% Aufgabe 10b
Beo_Matrix2 = [C; C*A; C*A*A; C*A*A*A];

rang_Beo_Matrix2 = rank(Beo_Matrix2);

%% Aufgabe 11

Ste_Matrix = [b, A*b, A*A*b, A*A*A*b];

det_Ste_Matrix = simplify(det(Ste_Matrix));

%% Aufgabe 13
A = subs(A, {g, m_K, m_L, l, J, R, K_E, K_M, K_U, K_red, r, m_K_schlange}, {9.81, 13.5, 5, 1, 0.13e-3, 0.5, 0.082, 0.08, 3.352, 21, 0.043, 44.5059});
b = subs(b, {g, m_K, m_L, l, J, R, K_E, K_M, K_U, K_red, r, m_K_schlange}, {9.81, 13.5, 5, 1, 0.13e-3, 0.5, 0.082, 0.08, 3.352, 21, 0.043, 44.5059});

A = double(vpa(A));
b = double(vpa(b));

Eigenwerte = eig(A);

%% Aufgabe 15

p = [1 2 3 4]; % Pole hinzufügen

k = acker(A,b,p); % Pegler

vorFilter = -inv(C_a*inv(A-b*k)*b);

%% Aufgabe 18

% Beobachter mit Place (Plovorgabe)
l = [1 2 3 4]; % Pole des Beobachters hinzufügen

L_place = (place(A',C',l))';

% Beobachter mit Riccati-Gleichung/ Kalman Filter

Q = eye(min(size(A)));
R = eye(min(size(C)));
E = eye(min(size(A)));
S = zeros(min(size(A)), min(size(C')));

[~,~,G_reccati] = care(A', C', Q, R, S, E);
L_riccati = G_reccati';

% L_voll = L_place oder L_riccati
L_voll = L_riccati;

%% Aufgabe 20
% Koordinatentransformation
T = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1];
Ar = inv(T)*A*T;
br = inv(T)*b;
Cr = C*T;
n = 2;

Ar11 = Ar(1:n,1:n);
Ar12 = Ar(1:n,n+1:length(A));
Ar21 = Ar(n+1:length(A), 1:n);
Ar22 = Ar(n+1:length(A),n+1:length(A));

br1 = br(1:n);
br2 = br(n+1:length(A));

Cr1 = Cr(1:n,1:n);
Cr2 = Cr(1:n, n+1:length(Cr));

% Beobachterentwurf mit Place/Polvorgabe
l = [1 2]; %Pole des Beobachters hinzufügen

L_place_red = (place(Ar22', Ar12', l))';

% Beobachterentwurf mit Place/Polvorgabe
Q_red = eye(min(size(Ar22)));
R_red = eye(min(size(Cr1)));
E_red = eye(min(size(Ar22)));
S_red = zeros(min(size(Ar22)), min(size(Cr1')));

[~,~,G_riccati_red] = care(Ar22', Cr1', Q_red, R_red, S_red, E_red);
L_riccati_red = G_riccati_red';

% L_red = L_place_red oder L_riccati_red;
L_red = L_riccati_red;










