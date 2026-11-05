%% RBF 3D
clear; clc; close all;

% --- Fonction cible 3D ---
f = @(x,y,z) sin(pi*x) .* cos(pi*y) .* exp(-z);

% --- Points de données dispersés en 3D ---
N = 200;
rng(42);
x_data = rand(N,1);
y_data = rand(N,1);
z_data = rand(N,1);
f_data = f(x_data, y_data, z_data);

% --- Paramètre de forme ---
epsilon = 3;

% --- Fonction RBF (Gaussienne) ---
phi = @(r) exp(-(epsilon*r).^2);

% --- Construction de la matrice A ---
A = zeros(N, N);
for i = 1:N
    for j = 1:N
        r = sqrt((x_data(i)-x_data(j))^2 + ...
                 (y_data(i)-y_data(j))^2 + ...
                 (z_data(i)-z_data(j))^2);
        A(i,j) = phi(r);
    end
end

% --- Résolution ---
w = A \ f_data;

% --- Évaluation sur une grille 3D (coupe à z=0.5) ---
[Xg, Yg] = meshgrid(linspace(0,1,50), linspace(0,1,50));
z_coupe = 0.5;
x_eval = Xg(:);
y_eval = Yg(:);
z_eval = z_coupe * ones(size(x_eval));

f_rbf = zeros(size(x_eval));
for i = 1:N
    r = sqrt((x_eval - x_data(i)).^2 + ...
             (y_eval - y_data(i)).^2 + ...
             (z_eval - z_data(i)).^2);
    f_rbf = f_rbf + w(i) * phi(r);
end
F_rbf = reshape(f_rbf, size(Xg));
F_exact = f(Xg, Yg, z_coupe);

% --- Affichage ---
figure;
subplot(1,2,1);
surf(Xg, Yg, F_exact); title('Fonction exacte (z=0.5)');
shading interp; colorbar; xlabel('x'); ylabel('y');

subplot(1,2,2);
surf(Xg, Yg, F_rbf); title('RBF 3D interpolée (z=0.5)');
shading interp; colorbar; xlabel('x'); ylabel('y');
sgtitle('RBF 3D - Coupe à z = 0.5');

% --- Erreur max ---
fprintf('Erreur max 3D : %.6f\n', max(abs(F_rbf(:) - F_exact(:))))