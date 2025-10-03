figure; 
t = 0:0.001:1;

decasteljau_check = uicontrol('Style', 'checkbox', ...
    'String', 'de Casteljau', ...
    'Value', 1, ...
    'Units', 'normalized', ...
    'Position', [0.8 0.02 0.2 0.05], ...
    'Callback', @updatePlot);

bernstein_check = uicontrol('Style', 'checkbox', ...
    'String', 'Bernstein', ...
    'Value', 0, ...
    'Units', 'normalized', ...
    'Position', [0.8 0.08 0.2 0.05], ...
    'Callback', @updatePlot);

h = drawpolyline;
hold on;

curve_points = bezier_decasteljau(h.Position, t);
hcurve1 = plot(curve_points(:,1), curve_points(:,2), 'g', 'linewidth', 2);
hcurve2 = plot(NaN, NaN, 'b', 'linewidth', 2);

h.addlistener('MovingROI', @(src, evt) updatePlot());
axis manual;

%%

function updatePlot(~, ~)
    control_points = evalin('base', 'h.Position');
    decasteljau_val = evalin('base', 'decasteljau_check.Value');
    bernstein_val = evalin('base', 'bernstein_check.Value');
    hcurve1 = evalin('base', 'hcurve1');
    hcurve2 = evalin('base', 'hcurve2');
    t = evalin('base', 't');
    
    if size(control_points, 1) < 2
        error('At least 2 points')
    end

    if decasteljau_val
        curve_points = bezier_decasteljau(control_points, t);
        set(hcurve1, 'XData', curve_points(:,1), 'YData', curve_points(:,2));
    else
        set(hcurve1, 'XData', NaN, 'YData', NaN);
    end

    if bernstein_val
        curve_points = bezier_bernstein(control_points, t);
        set(hcurve2, 'XData', curve_points(:,1), 'YData', curve_points(:,2));
    else
        set(hcurve2, 'XData', NaN, 'YData', NaN);
    end
end

%%

function p_curve = bezier_decasteljau(p, t)
    n = size(p, 1);
    m = length(t);
    p_curve = zeros(m, 2);
    
    for k = 1:m
        points = p;
        for r = 1:n-1
            for i = 1:n-r
                points(i,:) = (1 - t(k)) * points(i,:) + t(k) * points(i+1,:);
            end
        end
        p_curve(k,:) = points(1,:);
    end
end

function p_curve = bezier_bernstein(p, t)
    n = size(p, 1) - 1;
    m = length(t);
    p_curve = zeros(m, 2);
    
    for k = 1:m
        B = zeros(1, n+1);
        for i = 0:n
            B(i+1) = nchoosek(n, i) * (t(k)^i) * ((1 - t(k))^(n - i));
        end
        p_curve(k,:) = B * p;
    end
end
