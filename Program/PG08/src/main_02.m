function main_02()
    % 获取单位球的基础控制点 (双二次)
    s2 = sqrt(2)/2;
    W = [1,  s2, 1; s2, 0.5, s2; 1,  s2, 1];
    Px = [1, 1, 0; 1, 1, 0; 0, 0, 0];
    Py = [0, 1, 1; 0, 1, 1; 0, 0, 0];
    Pz = [0, 0, 0; 1, 1, 1; 1, 1, 1];

    % 变换为椭球: 3x^2 + 2y^2 + z^2 = 1
    % 半轴: a = 1/sqrt(3), b = 1/sqrt(2), c = 1
    Px = Px * (1/sqrt(3));
    Py = Py * (1/sqrt(2));
    Pz = Pz * 1; 

    % 升阶 (Degree Elevation): 从双二次(2x2) -> 双三次(3x3)
    % 我们需要在 u 和 v 方向分别进行升阶

    % 将由 (WX, WY, WZ, W) 组成的 4D 控制点进行升阶
    [CPx_3, CPy_3, CPz_3, W_3] = elevate_degree_2_to_3(Px, Py, Pz, W);

    % 绘制双三次曲面
    plot_rational_bezier_surface(CPx_3, CPy_3, CPz_3, W_3, [1 0 0]);
    title({'Problem 2:', 'Ellipsoid 3x^2+2y^2+z^2=1', 'Bicubic (Degree Elevated)'});
    axis equal; grid on; view(3);
    xlabel('X'); ylabel('Y'); zlabel('Z');
end

function [Nx, Ny, Nz, Nw] = elevate_degree_2_to_3(Ox, Oy, Oz, Ow)
    % 升阶函数：将 2次 Bezier 控制网格提升为 3次
    % 输入是 3x3 矩阵，输出是 4x4 矩阵

    % 组合成 4D 齐次坐标 (Homogeneous Coordinates): H = [wx, wy, wz, w]
    % 必须在齐次空间进行升阶！
    H = zeros(3, 3, 4);
    H(:,:,1) = Ox .* Ow; H(:,:,2) = Oy .* Ow; 
    H(:,:,3) = Oz .* Ow; H(:,:,4) = Ow;

    % 两次调用一维升阶：先对列(u方向)，再对行(v方向)
    H_temp = zeros(4, 3, 4); % u 方向升阶后: 3->4 行
    for col = 1:3
        for dim = 1:4
            H_temp(:, col, dim) = elevate_1d(H(:, col, dim));
        end
    end

    H_final = zeros(4, 4, 4); % v 方向升阶后: 3->4 列
    for row = 1:4
        for dim = 1:4
            H_final(row, :, dim) = elevate_1d(H_temp(row, :, dim)')';
        end
    end

    % 还原回欧几里得坐标
    Nw = H_final(:,:,4);
    Nx = H_final(:,:,1) ./ Nw;
    Ny = H_final(:,:,2) ./ Nw;
    Nz = H_final(:,:,3) ./ Nw;
end

function new_vec = elevate_1d(old_vec)
    % 1D Bézier 升阶: 2次 (3点) -> 3次 (4点)
    % 公式: Q_i = (i/n+1)*P_{i-1} + (1 - i/n+1)*P_i
    % n=2 (原次数), 新次数 3. i=0..3
    % old_vec = [P0, P1, P2]
    P = old_vec;
    Q = zeros(4, 1);
    Q(1) = P(1); % Q0 = P0
    Q(2) = (1/3)*P(1) + (2/3)*P(2);
    Q(3) = (2/3)*P(2) + (1/3)*P(3);
    Q(4) = P(3); % Q3 = P2
    new_vec = Q;
end
