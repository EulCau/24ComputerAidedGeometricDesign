function assignment7_q1()
    clc; close all;

    % ===== 用户输入 =====
    cx = input('请输入立方体中心 cx: ');
    cy = input('请输入立方体中心 cy: ');
    cz = input('请输入立方体中心 cz (>0): ');
    d  = input('请输入立方体半边长 d: ');

    % 焦距 (随便假设一个合理值)
    f = 1.0;

    % ===== 构造立方体顶点 (世界坐标 = 摄像机坐标) =====
    C = [cx; cy; cz];
    offsets = d * [
        -1 -1 -1;
        -1 -1  1;
        -1  1 -1;
        -1  1  1;
         1 -1 -1;
         1 -1  1;
         1  1 -1;
         1  1  1;
    ];
    V = offsets + C';     % 8x3 matrix, 每行一个顶点 [x y z]

    % ===== 透视投影 =====
    X = V(:,1);
    Y = V(:,2);
    Z = V(:,3);

    if any(Z <= 0)
        warning('有顶点 z <= 0, 可能会导致投影错误, 请调整中心 cz.');
    end

    U = f * X ./ Z;
    W = f * Y ./ Z;

    % 立方体的12条边（顶点索引）
    edges = [
        1 2; 1 3; 1 5;
        2 4; 2 6;
        3 4; 3 7;
        4 8;
        5 6; 5 7;
        6 8;
        7 8
    ];

    % ===== 画二维投影 =====
    figure; hold on; axis equal;
    for k = 1:size(edges,1)
        i = edges(k,1);
        j = edges(k,2);
        plot([U(i), U(j)], [W(i), W(j)], '-o');
    end
    xlabel('u'); ylabel('v');
    title('立方体的透视投影 (2D)');
    grid on;
end
