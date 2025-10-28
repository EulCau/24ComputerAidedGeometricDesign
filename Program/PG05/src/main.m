function main()
    close all;
    fig = figure('Name','Cubic B-Spline Editor','NumberTitle','off');
    ax = axes('Parent', fig);
    axis(ax, [0 1 0 1]);
    axis(ax, 'equal'); axis(ax, 'manual');
    grid(ax, 'on'); hold(ax, 'on');
    title(ax, {'点击添加锚点并回车结束','拖动锚点调整 B 样条曲线'});
    xlabel(ax, 'x'); ylabel(ax, 'y');

    disp('请点击添加锚点 (回车结束)...');
    roi = drawpolyline('Color','k','LineWidth',1.5);
    P = roi.Position;
    if size(P,1) < 4
        error('至少需要 4 个锚点!');
    end

    % 绘制 B 样条曲线
    [curve,poly] = plotBSpline(P);

    % 拖动锚点时更新曲线
    addlistener(roi, 'MovingROI', @(src, evt) updateBSpline(src, curve, poly));
end

%%

function [h,p] = plotBSpline(P)
    n = size(P,1);
    k = 4; % order (degree = 3)

    u = zeros(n,1);
    for i = 2:n
        u(i) = u(i-1) + norm(P(i,:) - P(i-1,:));
    end
    u = u / u(end);

    t = [zeros(1,k), 1:(n-k), (n-k+1)*ones(1,k)];
    t = t / (n-k+1);

    A = zeros(n);
    for j = 1:n
        for i = 1:n
            A(j,i) = bspline_basis(i, k, t, u(j));
        end
    end

    Q = A \ P;

    ts = linspace(0,1,300);
    B = zeros(length(ts),2);
    for ii = 1:length(ts)
        B(ii,:) = bspline_point(Q, k, t, ts(ii));
    end

    p = plot(Q(:,1), Q(:,2), 'ro--', 'LineWidth', 1.0);
    h = plot(B(:,1), B(:,2), 'b-', 'LineWidth', 2);
end


%%

function updateBSpline(roi, h, p)
    P = roi.Position;
    if size(P,1) < 4, return; end
    n = size(P,1);
    k = 4;

    u = zeros(n,1);
    for i = 2:n
        u(i) = u(i-1) + norm(P(i,:) - P(i-1,:));
    end
    u = u / u(end);

    t = [zeros(1,k), 1:(n-k), (n-k+1)*ones(1,k)];
    t = t / (n-k+1);

    A = zeros(n);
    for j = 1:n
        for i = 1:n
            A(j,i) = bspline_basis(i, k, t, u(j));
        end
    end
    Q = A \ P;

    ts = linspace(0,1,300);
    B = zeros(length(ts),2);
    for ii = 1:length(ts)
        B(ii,:) = bspline_point(Q, k, t, ts(ii));
    end

    set(h, 'XData', B(:,1), 'YData', B(:,2));
    set(p, 'XData', Q(:,1), 'YData', Q(:,2));
end

%%

function C = bspline_point(P, k, t, u)
    n = size(P,1);
    C = [0,0];
    for i = 1:n
        Ni = bspline_basis(i, k, t, u);
        C = C + Ni * P(i,:);
    end
end

%%

function N = bspline_basis(i, k, t, u)
    if k == 1
        if (t(i) <= u && u < t(i+1)) || (u == t(end) && t(i) <= u && u <= t(i+1))
            N = 1;
        else
            N = 0;
        end
    else
        denom1 = t(i+k-1) - t(i);
        denom2 = t(i+k) - t(i+1);

        term1 = 0;
        term2 = 0;
        if denom1 ~= 0
            term1 = (u - t(i)) / denom1 * bspline_basis(i, k-1, t, u);
        end
        if denom2 ~= 0
            term2 = (t(i+k) - u) / denom2 * bspline_basis(i+1, k-1, t, u);
        end
        N = term1 + term2;
    end
end
