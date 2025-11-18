function assignment7_q3()
    clc; close all;

    % 用户输入
    a = input('请输入椭圆/双曲线参数 a: ');
    b = input('请输入椭圆/双曲线参数 b: ');

    % ================= 椭圆在齐次坐标中的 3D 表示 =================
    figure; hold on; grid on;
    title('椭圆：齐次坐标中的3D Bézier曲线');
    xlabel('X (wx)'); ylabel('Y (wy)'); zlabel('W (权重)');
    view(3);

    % 正确的权重计算: 对于 90 度弧段, 权重 w = cos(45°) = √2/2
    w = sqrt(2)/2;

    % 四段椭圆弧的配置 (每段 90 度)
    segments = {
        % 第一象限 (0° to 90°)
        struct('P0', [a, 0], 'P1', [a, b], 'P2', [0, b], 'w', w);
        % 第二象限 (90° to 180°)  
        struct('P0', [0, b], 'P1', [-a, b], 'P2', [-a, 0], 'w', w);
        % 第三象限 (180° to 270°)
        struct('P0', [-a, 0], 'P1', [-a, -b], 'P2', [0, -b], 'w', w);
        % 第四象限 (270° to 360°)
        struct('P0', [0, -b], 'P1', [a, -b], 'P2', [a, 0], 'w', w);
    };

    for k = 1:length(segments)
        seg = segments{k};

        % 提取控制点和权重
        P0_2D = seg.P0';
        P1_2D = seg.P1';
        P2_2D = seg.P2';
        w1 = seg.w;
        w0 = 1;
        w2 = 1;

        % 将 2D 控制点转换为齐次坐标中的 3D 控制点
        % 在齐次坐标中, 表示为 (w*x, w*y, w)
        P0_3D = [w0 * P0_2D(1), w0 * P0_2D(2), w0];
        P1_3D = [w1 * P1_2D(1), w1 * P1_2D(2), w1];
        P2_3D = [w2 * P2_2D(1), w2 * P2_2D(2), w2];

        % 在齐次坐标中评估 3D Bézier 曲线 (非有理)
        t = linspace(0, 1, 100);
        B0 = (1-t).^2;
        B1 = 2*(1-t).*t;
        B2 = t.^2;

        % 3D Bézier曲线计算 (在齐次坐标空间中)
        x_3D = B0*P0_3D(1) + B1*P1_3D(1) + B2*P2_3D(1);
        y_3D = B0*P0_3D(2) + B1*P1_3D(2) + B2*P2_3D(2);
        w_3D = B0*P0_3D(3) + B1*P1_3D(3) + B2*P2_3D(3);

        % 绘制 3D 齐次坐标中的 Bézier 曲线
        plot3(x_3D, y_3D, w_3D, 'b-', 'LineWidth', 2);

        % 绘制控制点在齐次坐标中的位置
        plot3([P0_3D(1), P1_3D(1), P2_3D(1)], ...
              [P0_3D(2), P1_3D(2), P2_3D(2)], ...
              [P0_3D(3), P1_3D(3), P2_3D(3)], 'ro-', 'MarkerSize', 6);

        % 标注控制点
        text(P0_3D(1), P0_3D(2), P0_3D(3), 'P0', 'FontSize', 10);
        text(P1_3D(1), P1_3D(2), P1_3D(3), 'P1', 'FontSize', 10);
        text(P2_3D(1), P2_3D(2), P2_3D(3), 'P2', 'FontSize', 10);
    end

    % 添加投影平面 (w=1平面)
    [X, Y] = meshgrid(linspace(-a, a, 10), linspace(-b, b, 10));
    Z = ones(size(X));
    surf(X, Y, Z, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', 'green');

    legend('齐次坐标中的Bézier曲线', '控制多边形', '投影平面 (w=1)');

    % ================= 双曲线在齐次坐标中的 3D 表示 =================
    figure; hold on; grid on;
    title('双曲线右支：齐次坐标中的3D Bézier曲线');
    xlabel('X (wx)'); ylabel('Y (wy)'); zlabel('W (权重)');
    view(3);

    % 选择双曲线的一段进行表示
    t_start = 0.5;
    t_end = 2.0;

    % 计算双曲线上的点
    P0_2D = [a*cosh(t_start); b*sinh(t_start)];
    P2_2D = [a*cosh(t_end); b*sinh(t_end)];
    P1_mid_2D = [a*cosh((t_start+t_end)/2); b*sinh((t_start+t_end)/2)];

    % 计算权重
    theta = (t_end - t_start) / 2;
    w1 = cosh(theta);
    w0 = 1;
    w2 = 1;

    % 计算中间控制点
    t_mid = 0.5;
    B0_mid = (1-t_mid)^2;
    B1_mid = 2*(1-t_mid)*t_mid;
    B2_mid = t_mid^2;

    P1_2D = (P1_mid_2D * (B0_mid*w0 + B1_mid*w1 + B2_mid*w2) - ...
             B0_mid*w0*P0_2D - B2_mid*w2*P2_2D) / (B1_mid*w1);

    % 将 2D 控制点转换为齐次坐标中的3D控制点
    P0_3D = [w0 * P0_2D(1), w0 * P0_2D(2), w0];
    P1_3D = [w1 * P1_2D(1), w1 * P1_2D(2), w1];
    P2_3D = [w2 * P2_2D(1), w2 * P2_2D(2), w2];

    % 在齐次坐标中评估 3D Bézier 曲线
    t = linspace(0, 1, 100);
    B0 = (1-t).^2;
    B1 = 2*(1-t).*t;
    B2 = t.^2;

    x_3D = B0*P0_3D(1) + B1*P1_3D(1) + B2*P2_3D(1);
    y_3D = B0*P0_3D(2) + B1*P1_3D(2) + B2*P2_3D(2);
    w_3D = B0*P0_3D(3) + B1*P1_3D(3) + B2*P2_3D(3);

    % 绘制 3D 齐次坐标中的 Bézier 曲线
    plot3(x_3D, y_3D, w_3D, 'b-', 'LineWidth', 2);

    % 绘制控制点在齐次坐标中的位置
    plot3([P0_3D(1), P1_3D(1), P2_3D(1)], ...
          [P0_3D(2), P1_3D(2), P2_3D(2)], ...
          [P0_3D(3), P1_3D(3), P2_3D(3)], 'ro-', 'MarkerSize', 6);

    % 标注控制点
    text(P0_3D(1), P0_3D(2), P0_3D(3), 'P0', 'FontSize', 10);
    text(P1_3D(1), P1_3D(2), P1_3D(3), 'P1', 'FontSize', 10);
    text(P2_3D(1), P2_3D(2), P2_3D(3), 'P2', 'FontSize', 10);

    % 添加投影平面 (w=1 平面)
    x_range = linspace(P0_2D(1), P2_2D(1), 10);
    y_range = linspace(min(P0_2D(2), P2_2D(2)), max(P0_2D(2), P2_2D(2)), 10);
    [X, Y] = meshgrid(x_range, y_range);
    Z = ones(size(X));
    surf(X, Y, Z, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', 'green');

    legend('齐次坐标中的Bézier曲线', '控制多边形', '投影平面 (w=1)');

    % ================= 显示投影变换 =================
    % 创建一个新图形, 显示从 3D 齐次坐标到 2D 标准坐标的投影
    figure; 

    % 子图 1: 3D 齐次坐标
    subplot(1,2,1); hold on; grid on;
    title('3D齐次坐标空间');
    xlabel('wx'); ylabel('wy'); zlabel('w');
    view(3);

    % 绘制椭圆的一段在 3D 齐次坐标中
    seg = segments{1}; % 使用第一段椭圆
    P0_2D = seg.P0';
    P1_2D = seg.P1';
    P2_2D = seg.P2';
    w1 = seg.w;
    w0 = 1;
    w2 = 1;

    P0_3D = [w0 * P0_2D(1), w0 * P0_2D(2), w0];
    P1_3D = [w1 * P1_2D(1), w1 * P1_2D(2), w1];
    P2_3D = [w2 * P2_2D(1), w2 * P2_2D(2), w2];

    t = linspace(0, 1, 100);
    B0 = (1-t).^2;
    B1 = 2*(1-t).*t;
    B2 = t.^2;

    x_3D = B0*P0_3D(1) + B1*P1_3D(1) + B2*P2_3D(1);
    y_3D = B0*P0_3D(2) + B1*P1_3D(2) + B2*P2_3D(2);
    w_3D = B0*P0_3D(3) + B1*P1_3D(3) + B2*P2_3D(3);

    plot3(x_3D, y_3D, w_3D, 'b-', 'LineWidth', 2);
    plot3([P0_3D(1), P1_3D(1), P2_3D(1)], ...
          [P0_3D(2), P1_3D(2), P2_3D(2)], ...
          [P0_3D(3), P1_3D(3), P2_3D(3)], 'ro-', 'MarkerSize', 6);

    % 添加投影线 (从 3D 点到投影平面)
    for i = 1:10:length(x_3D)
        plot3([x_3D(i), x_3D(i)/w_3D(i)], ...
              [y_3D(i), y_3D(i)/w_3D(i)], ...
              [w_3D(i), 1], 'k--', 'LineWidth', 0.5);
    end

    % 子图 2: 2D 投影结果
    subplot(1,2,2); hold on; axis equal; grid on;
    title('2D投影结果 (除以w坐标)');
    xlabel('x'); ylabel('y');

    % 计算 2D 投影点
    x_2D = x_3D ./ w_3D;
    y_2D = y_3D ./ w_3D;

    plot(x_2D, y_2D, 'b-', 'LineWidth', 2);

    % 绘制 2D 控制点
    P0_2D_proj = P0_3D(1:2) / P0_3D(3);
    P1_2D_proj = P1_3D(1:2) / P1_3D(3);
    P2_2D_proj = P2_3D(1:2) / P2_3D(3);

    plot([P0_2D_proj(1), P1_2D_proj(1), P2_2D_proj(1)], ...
         [P0_2D_proj(2), P1_2D_proj(2), P2_2D_proj(2)], 'ro-', 'MarkerSize', 6);

    fprintf('第三题完成:\n');
    fprintf('1. 显示了椭圆和双曲线在齐次坐标中的3D Bézier曲线\n');
    fprintf('2. 控制点在3D空间中的坐标为 (w*x, w*y, w)\n');
    fprintf('3. 通过除以w坐标（投影到w=1平面）可以得到2D标准坐标中的曲线\n');
end
