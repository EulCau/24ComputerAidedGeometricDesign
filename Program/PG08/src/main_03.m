function main_03()
    % 控制点
    P = [0 0 0; 2 2 4; 4 -2 6; 4 4 0; 8 0 4; 6 -4 4]'; % 3x6

    % 注意：三角形网格连线比较复杂，这里只画控制点
    plot3(P(1,:), P(2,:), P(3,:), 'rs', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    text(P(1,1), P(2,1), P(3,1), ' a(00)');
    hold on;

    % 待测参数
    uv_list = [0.25 0.5; 0.3 0.75; 0.5 0.5];
    names = {'p1', 'p2', 'p3'};

    for k = 1:3
        u = uv_list(k, 1);
        v = uv_list(k, 2);

        % 计算重心坐标 (基于题意: c=(0.5,1), b=(1,0), a=(0,0))
        % explicit formulas derived previously:
        gamma = v;
        beta = u - 0.5*v;
        alpha = 1 - beta - gamma;

        fprintf('Point %s: (alpha=%.2f, beta=%.2f, gamma=%.2f) -> ', names{k}, alpha, beta, gamma);

        if alpha < -1e-6 || beta < -1e-6 || gamma < -1e-6
            fprintf('Outside Triangle\n');
        else
            fprintf('Inside Triangle\n');
            % De Casteljau
            pt = de_casteljau_tri(P, alpha, beta, gamma);
            plot3(pt(1), pt(2), pt(3), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
            text(pt(1), pt(2), pt(3), [' ' names{k}]);
            fprintf('   Surface Point: (%.4f, %.4f, %.4f)\n', pt);
        end
    end

    title({'Problem 3:', 'Bezier Triangle', '& Point Evaluation'});
    axis equal; grid on; view(3);
    hold off;
end

function pt = de_casteljau_tri(P, a, b, c)
    % r=1
    P1_100 = a*P(:,1) + b*P(:,2) + c*P(:,3);
    P1_010 = a*P(:,2) + b*P(:,4) + c*P(:,5);
    P1_001 = a*P(:,3) + b*P(:,5) + c*P(:,6);
    % r=2
    pt = a*P1_100 + b*P1_010 + c*P1_001;
end
