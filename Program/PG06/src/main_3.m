%% Assignment 6 - Problem 3 (Robust de Boor + Bezier conversion)
% ----------------------------------------------------------
% This script evaluates a uniform cubic B-spline at t=2.5 using
% a robust, span-consistent de Boor implementation, and converts
% the spline segment [2,3] into an equivalent cubic Bezier curve.
% All comments are in English as requested.
% ----------------------------------------------------------

clear; clc; close all;

% --- Given control points and knot vector (uniform, open) ---
P = [-2 -10; -4  2;  6  5;  4 -7];
U = [ 0   0   1   2   3   4   5   5];
p = 3;                  % degree (cubic)
t_eval = 2.5;           % evaluation parameter (inside [2,3])

%% ================= Helper: findKnotSpan (MATLAB 1-based) =================
% Returns k such that U(k) <= t < U(k+1), with k in [p+1, n+1],
% where n = number_of_ctrl_pts - 1.
function k = findKnotSpan(U, p, t)
    m = numel(U) - 1;           % last knot index (0-based notion -> here "m")
    n = m - p - 1;              % n = #ctrlPts - 1
    % Right-end protection: if t == last internal knot U(n+2), put k = n+1
    if abs(t - U(n+2)) < 1e-14
        k = n + 1;
        return;
    end
    % Linear search is fine for small vectors
    for kk = (p+1):(n+1)
        if t >= U(kk) && t < U(kk+1)
            k = kk;
            return;
        end
    end
    error('t is out of the valid parametric range for this spline segment.');
end

%% ================= Helper: de Boor (robust, span-consistent) ==============
% Returns the curve point C at parameter t, and 'layers' for visualization.
% The algorithm uses a local (p+1)-by-d array d, consistent with the span k.
function [C, layers] = deBoor_full(U, P, p, t)
    % Determine sizes
    d_dim = size(P,2);                 % dimension (2D here)
    % Find span index k such that U(k) <= t < U(k+1)
    k = findKnotSpan(U, p, t);
    % Initialize local working array of length p+1
    d = zeros(p+1, d_dim);
    % Copy relevant control points: P(k-p : k)
    % NOTE: MATLAB is 1-based; indices k-p ... k are valid for k in [p+1, n+1]
    for j = 0:p
        d(j+1,:) = P(k - p + j, :);
    end
    % Store the first "layer" (optional: for plotting the construction)
    layers = cell(p+1,1);
    layers{1} = d;  % these are the p+1 control points involved

    % De Boor recursion
    for r = 1:p
        for j = p:-1:r
            denom = U(k + 1 + j - r) - U(k - p + j);
            % Numerical safety: if denom is zero (shouldn't happen for valid span),
            % clamp alpha to 0 (or 1) to avoid NaN. Here we just guard:
            if abs(denom) < 1e-14
                alpha = 0.0;
            else
                alpha = (t - U(k - p + j)) / denom;
            end
            d(j+1,:) = (1 - alpha)*d(j,:) + alpha*d(j+1,:);
        end
        % Keep only the "active" part for this layer (first p-r+1 points)
        layers{r+1} = d(1:(p - r + 1), :);
    end
    C = d(p+1,:);   % final point
end

%% ================= (1) Evaluate at t = 2.5 via de Boor ====================
[C, layers] = deBoor_full(U, P, p, t_eval);
disp('Curve point at t = 2.5:');
disp(C);   % Expected: [1.0000, 3.0000]

% --- Plot de Boor construction ---
figure; hold on; axis equal; grid on;
title('Problem 3: de Boor construction at t = 2.5');
xlabel('x'); ylabel('y');

% Original control polygon
plot(P(:,1), P(:,2), 'ko--', 'LineWidth', 1.2, ...
     'DisplayName', 'B-spline control polygon');

% Plot intermediate layers (local points in each recursion)
colors = {'r','g','b','m'};
for r = 2:length(layers)
    pts = layers{r};
    plot(pts(:,1), pts(:,2), 'o-', 'Color', colors{r-1}, ...
         'DisplayName', sprintf('Layer r=%d', r-1));
end

% Final curve point
plot(C(1), C(2), 'ks', 'MarkerFaceColor', 'k', ...
     'DisplayName', 'C(2.5)');

legend('Location','best');

%% ================= (2) Convert [2,3] to cubic Bezier ======================
% We need C(2), C'(2), C(3), C'(3). We'll evaluate C at the endpoints,
% and approximate derivatives numerically with a small step h. To avoid
% evaluating exactly at the right endpoint knot (which may produce 0 denominators),
% we use "3 - h" instead of "3" directly.

h = 1e-6;

% Endpoints (avoid exact right-end knot)
C2   = deBoor_full(U, P, p, 2.0);
C3   = deBoor_full(U, P, p, 3.0 - h);

% First derivatives via one-sided finite differences
C2p  = (deBoor_full(U, P, p, 2.0 + h) - C2) / h;
C3p  = (C3 - deBoor_full(U, P, p, 3.0 - 2*h)) / h;

% Bezier control points over length Delta = 1
B0 = C2;
B1 = C2 + (1/3) * C2p;
B2 = C3 - (1/3) * C3p;
B3 = C3;
B  = [B0; B1; B2; B3];

disp('Equivalent cubic Bezier control points [x y]:');
disp(B);

% --- Helper: de Casteljau for Bezier evaluation ---
function Q = deCasteljau(Pb, t)
    n = size(Pb,1);
    for rr = 1:n-1
        for ii = 1:n-rr
            Pb(ii,:) = (1-t)*Pb(ii,:) + t*Pb(ii+1,:);
        end
    end
    Q = Pb(1,:);
end

% Plot Bezier control polygon and curve (should match the spline segment)
tt = linspace(0,1,200);
curve = zeros(numel(tt),2);
for i = 1:numel(tt)
    curve(i,:) = deCasteljau(B, tt(i));
end

plot(B(:,1), B(:,2), 'c--o', 'LineWidth', 1.2, ...
     'DisplayName', 'Bezier control polygon');
plot(curve(:,1), curve(:,2), 'c', 'LineWidth', 1.8, ...
     'DisplayName', 'Bezier curve');

legend('Location','best');
