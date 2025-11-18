function assignment7_q2()
    clc; close all;

    % 用户输入
    a = input('请输入椭圆/双曲线参数 a: ');
    b = input('请输入椭圆/双曲线参数 b: ');

    % ================= 椭圆 (有理二次 Bézier) =================
    figure; hold on; axis equal;
    title('椭圆：有理二次 Bézier 段拼接');
    xlabel('x'); ylabel('y');
    grid on;

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
        P0 = seg.P0';
        P1 = seg.P1';
        P2 = seg.P2';
        w1 = seg.w;
        w0 = 1;
        w2 = 1;

        % 评估有理二次 Bézier 曲线
        t = linspace(0, 1, 100);
        B0 = (1-t).^2;
        B1 = 2*(1-t).*t;
        B2 = t.^2;

        % 有理 Bézier 计算
        x_tilde = B0*w0*P0(1) + B1*w1*P1(1) + B2*w2*P2(1);
        y_tilde = B0*w0*P0(2) + B1*w1*P1(2) + B2*w2*P2(2);
        w_tilde = B0*w0 + B1*w1 + B2*w2;

        x = x_tilde ./ w_tilde;
        y = y_tilde ./ w_tilde;

        h1 = plot(x, y, 'b-', 'LineWidth', 2);
    end

    % 画理论椭圆作为对比
    t0 = linspace(0, 2*pi, 400);
    xe = a*cos(t0);
    ye = b*sin(t0);
    h2 = plot(xe, ye, 'r--', 'LineWidth', 1);
    legend([h1, h2], {'有理二次 Bézier 椭圆','理论椭圆'});

    % ================= 双曲线 (正确的有理二次 Bézier 表示) =================
    figure; hold on; axis equal;
    title('双曲线右支: 有理二次 Bézier 表示');
    xlabel('x'); ylabel('y');
    grid on;

    % 对于双曲线 x^2/a^2 - y^2/b^2 = 1, 我们可以用有理二次Bézier精确表示无穷大的一段
    % 但通常我们只表示有限的一段

    % 选择要绘制的双曲线范围 (参数 t 的范围)
    t_start = 0.5;
    t_end = 2.0;

    % 使用双曲线的标准参数方程: x = a*cosh(t), y = b*sinh(t)
    t_param = linspace(t_start, t_end, 3); % 只需要三个点来定义二次Bézier

    % 端点
    P0 = [a*cosh(t_param(1)); b*sinh(t_param(1))];
    P2 = [a*cosh(t_param(3)); b*sinh(t_param(3))];

    % 中间点 (在双曲线上)
    P1_mid = [a*cosh(t_param(2)); b*sinh(t_param(2))];

    % 对于双曲线的有理二次 Bézier 表示, 我们需要计算正确的权重
    % 根据圆锥曲线的有理表示理论, 权重与角度有关
    theta = (t_param(3) - t_param(1)) / 2;
    w1 = cosh(theta);  % 双曲线使用cosh而不是cos

    w0 = 1;
    w2 = 1;

    % 控制点 P1 需要通过权重调整位置
    % 在齐次坐标中, P1 = (w1 * P1_homogeneous) / w1
    % 所以我们需要找到合适的 P1_homogeneous 使得曲线通过 P1_mid

    % 在 t=0.5 时, 曲线应该通过 P1_mid
    t_mid = 0.5;
    B0_mid = (1-t_mid)^2;
    B1_mid = 2*(1-t_mid)*t_mid;
    B2_mid = t_mid^2;

    % 解出 P1 的位置
    P1 = (P1_mid * (B0_mid*w0 + B1_mid*w1 + B2_mid*w2) - ...
          B0_mid*w0*P0 - B2_mid*w2*P2) / (B1_mid*w1);

    % 评估有理二次 Bézier 曲线
    t = linspace(0, 1, 100);
    B0 = (1-t).^2;
    B1 = 2*(1-t).*t;
    B2 = t.^2;

    x_tilde = B0*w0*P0(1) + B1*w1*P1(1) + B2*w2*P2(1);
    y_tilde = B0*w0*P0(2) + B1*w1*P1(2) + B2*w2*P2(2);
    w_tilde = B0*w0 + B1*w1 + B2*w2;

    x_bezier = x_tilde ./ w_tilde;
    y_bezier = y_tilde ./ w_tilde;

    h1 = plot(x_bezier, y_bezier, 'b-', 'LineWidth', 2);

    % 理论双曲线 (相同范围)
    t_h = linspace(t_start, t_end, 200);
    xh = a * cosh(t_h);
    yh = b * sinh(t_h);
    h2 = plot(xh, yh, 'r--', 'LineWidth', 1);

    legend([h1, h2], {'有理二次 Bézier','理论双曲线'});
end
