%% RBF 3D
clear; clc; close all;
% --- 3D target function ---
f = @(x,y,z) sin(pi*x) .* cos(pi*y) .* exp(-z);
% --- Scattered data points in 3D ---
N = 200;
rng(42);
x_data = rand(N,1);
y_data = rand(N,1);
z_data = rand(N,1);
f_data = f(x_data, y_data, z_data);
% --- Shape parameter ---
epsilon = 3;
% --- RBF function (Gaussian) ---
phi = @(r) exp(-(epsilon*r).^2);
% --- Building matrix A ---
A = zeros(N, N);
for i = 1:N
for j = 1:N
        r = sqrt((x_data(i)-x_data(j))^2 + ...
                 (y_data(i)-y_data(j))^2 + ...
                 (z_data(i)-z_data(j))^2);
        A(i,j) = phi(r);
end
end
% --- Solving ---
w = A \ f_data;
% --- Evaluation on a 3D grid (slice at z=0.5) ---
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
% --- Display ---
figure;
subplot(1,2,1);
surf(Xg, Yg, F_exact); title('Exact function (z=0.5)');
shading interp; colorbar; xlabel('x'); ylabel('y');
subplot(1,2,2);
surf(Xg, Yg, F_rbf); title('3D interpolated RBF (z=0.5)');
shading interp; colorbar; xlabel('x'); ylabel('y');
sgtitle('3D RBF - Slice at z = 0.5');
% --- Max error ---
fprintf('Max error 3D: %.6f\n', max(abs(F_rbf(:) - F_exact(:))))
