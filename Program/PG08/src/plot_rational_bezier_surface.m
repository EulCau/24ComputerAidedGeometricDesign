function plot_rational_bezier_surface(Px, Py, Pz, W, color_code)
    [n_rows, n_cols] = size(W);
    n = n_rows - 1; % degree u
    m = n_cols - 1; % degree v

    u = linspace(0, 1, 30);
    v = linspace(0, 1, 30);
    [U, V] = meshgrid(u, v);

    SurfX = zeros(size(U)); SurfY = zeros(size(U)); SurfZ = zeros(size(U));

    % 计算曲面点 (de Casteljau 或 Bernstein polynomials)
    for i = 1:numel(U)
        uu = U(i); vv = V(i);

        % Bernstein 基函数
        Bu = bernstein_basis(n, uu);
        Bv = bernstein_basis(m, vv);

        num_x = 0; num_y = 0; num_z = 0; den = 0;

        for r = 0:n
            for c = 0:m
                w_val = W(r+1, c+1);
                b_val = Bu(r+1) * Bv(c+1);

                num_x = num_x + b_val * Px(r+1, c+1) * w_val;
                num_y = num_y + b_val * Py(r+1, c+1) * w_val;
                num_z = num_z + b_val * Pz(r+1, c+1) * w_val;
                den   = den   + b_val * w_val;
            end
        end

        SurfX(i) = num_x / den;
        SurfY(i) = num_y / den;
        SurfZ(i) = num_z / den;
    end

    surf(SurfX, SurfY, SurfZ, 'FaceAlpha', 0.7, 'FaceColor', color_code);
    hold on;
    % 绘制控制网格
    mesh(Px, Py, Pz, 'EdgeColor', 'k', 'FaceAlpha', 0);
    plot3(Px(:), Py(:), Pz(:), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
end

function B = bernstein_basis(n, t)
    B = zeros(n+1, 1);
    for i = 0:n
        B(i+1) = nchoosek(n, i) * (t^i) * ((1-t)^(n-i));
    end
end
