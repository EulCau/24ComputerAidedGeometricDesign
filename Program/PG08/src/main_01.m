function main_01()
    % 权重矩阵
    s2 = sqrt(2)/2;
    W = [1,  s2, 1;
         s2, 0.5, s2;
         1,  s2, 1];

    % X 坐标控制点
    Px = [1, 1, 0;
          1, 1, 0;
          0, 0, 0];

    % Y 坐标控制点
    Py = [0, 1, 1;
          0, 1, 1;
          0, 0, 0];

    % Z 坐标控制点
    Pz = [0, 0, 0;
          1, 1, 1;
          1, 1, 1];

    % 绘制曲面
    plot_rational_bezier_surface(Px, Py, Pz, W, [0 0 1]);
    title({'Problem 1:', 'Unit Sphere (Octant)', 'Biquadratic Rational'});
    axis equal; grid on; view(3);
    xlabel('X'); ylabel('Y'); zlabel('Z');
end
