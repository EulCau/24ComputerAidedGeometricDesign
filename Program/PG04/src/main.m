function main()
    close all;
    fig = figure('Name','Bézier Spline Editor','NumberTitle','off');
    ax = axes('Parent', fig);
    axis(ax, [0 1 0 1]);
    axis(ax, 'equal');
    axis(ax, 'manual');
    grid(ax, 'on'); hold(ax, 'on');
    title(ax, {'点击添加锚点并回车结束','拖动锚点或红点调整曲线','按 S 键切换 Smooth / Free 模式'});
    xlabel(ax, 'x'); ylabel(ax, 'y');

    % 记录当前模式
    mode = "Smooth";  % 默认光滑模式
    set(fig, 'UserData', mode);

    % 监听键盘切换
    set(fig, 'KeyPressFcn', @toggleMode);

    % 用户交互添加锚点
    disp('请点击添加锚点 (回车结束)...');
    roi = drawpolyline('Color','k','LineWidth',1.5);
    P = roi.Position;
    if size(P,1) < 2
        error('至少需要 2 个锚点!');
    end

    % 计算控制点
    [C1, C2] = computeBezierControls(P);
    n = size(C1,1);

    % 绘制控制点与曲线
    red1 = gobjects(n,1);
    red2 = gobjects(n,1);
    curve = gobjects(n,1);
    for i = 1:n
        red1(i) = drawpoint('Position', C1(i,:), 'Color', 'r');
        red2(i) = drawpoint('Position', C2(i,:), 'Color', 'r');
        curve(i) = plotBezier(P(i,:), C1(i,:), C2(i,:), P(i+1,:));
    end

    % 监听锚点拖动
    addlistener(roi, 'MovingROI', @(src, evt) updateAll(src, red1, red2, curve));

    % 监听控制点拖动
    for i = 1:n
        addlistener(red1(i), 'MovingROI', @(src,evt) updateCurve(fig, roi, red1, red2, curve, i, 1));
        addlistener(red2(i), 'MovingROI', @(src,evt) updateCurve(fig, roi, red1, red2, curve, i, 2));
    end
end

%% === Bézier 绘制函数 ===
function h = plotBezier(P0,P1,P2,P3)
    ts = linspace(0,1,100);
    B = (1-ts)'.^3*P0 + 3*(1-ts)'.^2.*ts'*P1 + ...
        3*(1-ts)'.*ts'.^2*P2 + ts'.^3*P3;
    h = plot(B(:,1), B(:,2), 'b-', 'LineWidth', 2);
end

%% === 控制点计算 ===
function [C1,C2] = computeBezierControls(P)
    n = size(P,1)-1;
    if n == 1
        C1 = (2*P(1,:)+P(2,:))/3;
        C2 = (P(1,:)+2*P(2,:))/3;
        return;
    end
    rhs = zeros(n,2);
    rhs(1,:) = P(1,:) + 2*P(2,:);
    for i=2:n-1
        rhs(i,:) = 4*P(i,:) + 2*P(i+1,:);
    end
    rhs(n,:) = 8*P(n,:) + P(n+1,:);
    b = [2;4*ones(n-2,1);7];
    C1 = zeros(n,2);
    for dim=1:2
        f = rhs(:,dim);
        bp=b;
        for i=2:n
            m=1/bp(i-1);
            bp(i)=bp(i)-m;
            f(i)=f(i)-m*f(i-1);
        end
        C1(n,dim)=f(n)/bp(n);
        for i=n-1:-1:1
            C1(i,dim)=(f(i)-C1(i+1,dim))/bp(i);
        end
    end
    C2=zeros(n,2);
    for i=1:n-1
        C2(i,:)=2*P(i+1,:)-C1(i+1,:);
    end
    C2(n,:)=(P(n+1,:)+C1(n,:))/2;
end

%% === 更新整条曲线（拖锚点） ===
function updateAll(roi, red1, red2, curve)
    P = roi.Position;
    n = numel(red1);
    for i=1:n
        p1 = red1(i).Position;
        p2 = red2(i).Position;
        updateBezier(curve(i), P(i,:), p1, p2, P(i+1,:));
    end
end

%% === 更新局部 Bézier 段 (拖控制点) ===
function updateCurve(fig, roi, red1, red2, curve, seg, which)
    mode = get(fig, 'UserData');
    P = roi.Position;
    n = numel(red1);

    % 当前段的控制点
    P1 = red1(seg).Position;
    P2 = red2(seg).Position;

    % 如果是 Smooth 模式: 保持共线性
    if mode == "Smooth"
        if which == 1 && seg > 1
            % 左控制点动, 更新上段的右控制点
            P_center = P(seg,:);
            new_dir = P_center - P1;      % 向量方向
            len = norm(red2(seg-1).Position - P(seg,:)); % 对侧长度
            if len > 0
                red2(seg-1).Position = P_center + len * new_dir / norm(new_dir);
            end
        elseif which == 2 && seg < n
            % 右控制点动, 更新下段的左控制点
            P_center = P(seg+1,:);
            new_dir = P_center - P2;
            len = norm(red1(seg+1).Position - P_center);
            if len > 0
                red1(seg+1).Position = P_center + len * new_dir / norm(new_dir);
            end
        end
    end

    % 更新整个曲线
    for i=1:n
        updateBezier(curve(i), P(i,:), red1(i).Position, red2(i).Position, P(i+1,:));
    end
end

%% === 仅更新曲线数据 ===
function updateBezier(h, P0,P1,P2,P3)
    ts = linspace(0,1,100);
    B = (1-ts)'.^3*P0 + 3*(1-ts)'.^2.*ts'*P1 + ...
        3*(1-ts)'.*ts'.^2*P2 + ts'.^3*P3;
    set(h, 'XData', B(:,1), 'YData', B(:,2));
end

%% === 模式切换函数 ===
function toggleMode(src, event)
    if event.Key == 's'
        current = get(src, 'UserData');
        if current == "Smooth"
            newMode = "Free";
        else
            newMode = "Smooth";
        end
        set(src, 'UserData', newMode);
        title(gca, sprintf('当前模式: %s (按 S 切换)', newMode));
        disp(['模式切换为: ', newMode]);
    end
end
