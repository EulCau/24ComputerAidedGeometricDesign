[x, t] = readObj('ball.obj');

figure; subplot(121); trimesh(t, x(:,1), x(:,2), x(:,3), 'edgecolor', 'k'); axis equal; axis off; title('coarse mesh');

[x2, t2] = loop(x, t, 1);

subplot(122); trimesh(t2, x2(:,1), x2(:,2), x2(:,3), 'edgecolor', 'k'); axis equal; axis off; title('subdivision mesh');


function [x2, t2] = loop(x, t, n)
    x2 = x; 
    t2 = t;

    for it = 1:n
        x_old = x2;
        t_old = t2;

        V = size(x_old, 1);   % 当前顶点数
        F = size(t_old, 1);   % 当前三角形数

        %% 1. 构建所有边及其“对顶点”信息
        % 每个三角形 (a,b,c) 产生三条边：
        % (a,b) 对顶点 c
        % (b,c) 对顶点 a
        % (c,a) 对顶点 b
        E_all = [t_old(:,[1 2]); t_old(:,[2 3]); t_old(:,[3 1])];   % (3F x 2)
        Opp   = [t_old(:,3);      t_old(:,1);      t_old(:,2)  ];   % (3F x 1)

        % 无向边：对每条边的两个端点排序，方便去重
        edges_sorted = sort(E_all, 2);

        % uniqueEdges: 所有不重复的无向边
        % idxEdge: E_all 每一行对应的 uniqueEdges 的索引
        [uniqueEdges, ~, idxEdge] = unique(edges_sorted, 'rows');
        E = size(uniqueEdges, 1);   % 总边数

        % 每条无向边最多有两个对顶点（内部边），边界边只有一个
        opp1 = zeros(E, 1);
        opp2 = zeros(E, 1);
        for r = 1:3*F
            e = idxEdge(r);
            if opp1(e) == 0
                opp1(e) = Opp(r);
            elseif opp2(e) == 0 && Opp(r) ~= opp1(e)
                opp2(e) = Opp(r);
            end
        end

        %% 2. 识别边界边 / 边界点，并收集边界邻居
        boundaryV = false(V, 1);       % 是否为边界点
        bNeigh    = cell(V, 1);        % 每个点的边界邻居（沿边界的两个点）

        for e = 1:E
            if opp2(e) == 0   % 只有一个对顶点 -> 边界边
                v0 = uniqueEdges(e,1);
                v1 = uniqueEdges(e,2);
                boundaryV(v0) = true;
                boundaryV(v1) = true;
                bNeigh{v0} = [bNeigh{v0}, v1];
                bNeigh{v1} = [bNeigh{v1}, v0];
            end
        end
        for i = 1:V
            if ~isempty(bNeigh{i})
                bNeigh{i} = unique(bNeigh{i});
            end
        end

        %% 3. 构建每个顶点的一环邻居（所有相邻点）
        neigh = cell(V, 1);
        for f = 1:F
            a = t_old(f,1); 
            b = t_old(f,2); 
            c = t_old(f,3);
            neigh{a} = [neigh{a}, b, c];
            neigh{b} = [neigh{b}, a, c];
            neigh{c} = [neigh{c}, a, b];
        end
        for i = 1:V
            if ~isempty(neigh{i})
                neigh{i} = unique(neigh{i});
            end
        end

        %% 4. 更新旧顶点位置（even vertices）
        x_even = zeros(V, 3);
        for i = 1:V
            if boundaryV(i)
                % 边界点规则：
                % v' = 3/4 * v + 1/8 * (v1 + v2)
                % 其中 v1, v2 为边界上的相邻两个点
                bn = bNeigh{i};
                if numel(bn) >= 2
                    v1 = bn(1);
                    v2 = bn(end);   % 一般就是两个
                    x_even(i,:) = 0.75 * x_old(i,:) + 0.125 * (x_old(v1,:) + x_old(v2,:));
                else
                    % 极端退化情况：就按内部点处理
                    Ni = neigh{i};
                    m = numel(Ni);
                    if m > 0
                        beta = (1/m) * (5/8 - ( (3/8) + (1/4)*cos(2*pi/m) )^2);
                        x_even(i,:) = (1 - m*beta)*x_old(i,:) + beta * sum(x_old(Ni,:), 1);
                    else
                        x_even(i,:) = x_old(i,:);
                    end
                end
            else
                % 内部点规则：
                % v' = (1 - n*β) * v + β * sum(v_i)
                % β = 1/n * [ 5/8 - ( 3/8 + 1/4 cos(2π/n) )^2 ]
                Ni = neigh{i};
                m = numel(Ni);
                if m > 0
                    beta = (1/m) * (5/8 - ( (3/8) + (1/4)*cos(2*pi/m) )^2);
                    x_even(i,:) = (1 - m*beta)*x_old(i,:) + beta * sum(x_old(Ni,:), 1);
                else
                    x_even(i,:) = x_old(i,:);
                end
            end
        end

        %% 5. 生成边上的新顶点（odd vertices）
        % 内部边：v = 3/8*(v0+v1) + 1/8*(v2+v3)
        % 边界边：v = 1/2*(v0+v1)
        x_odd = zeros(E, 3);
        for e = 1:E
            v0 = uniqueEdges(e,1);
            v1 = uniqueEdges(e,2);
            k1 = opp1(e);
            k2 = opp2(e);

            if k2 ~= 0
                % 内部边
                x_odd(e,:) = 3/8*(x_old(v0,:) + x_old(v1,:)) + ...
                             1/8*(x_old(k1,:) + x_old(k2,:));
            else
                % 边界边
                x_odd(e,:) = 0.5*(x_old(v0,:) + x_old(v1,:));
            end
        end

        %% 6. 组合新顶点数组
        % 先是更新后的老顶点，再接新生成的边顶点
        x_new = [x_even; x_odd];

        %% 7. 每个旧三角形拆成 4 个小三角形
        t_new = zeros(4*F, 3);
        for f = 1:F
            a = t_old(f,1);
            b = t_old(f,2);
            c = t_old(f,3);

            e1 = idxEdge(f);        % 对应 (a,b)
            e2 = idxEdge(F + f);    % 对应 (b,c)
            e3 = idxEdge(2*F + f);  % 对应 (c,a)

            ab = V + e1;
            bc = V + e2;
            ca = V + e3;

            base = 4*(f-1);
            t_new(base+1,:) = [a  ab ca];
            t_new(base+2,:) = [b  bc ab];
            t_new(base+3,:) = [c  ca bc];
            t_new(base+4,:) = [ab bc ca];
        end

        % 更新到下一层
        x2 = x_new;
        t2 = t_new;
    end
end
