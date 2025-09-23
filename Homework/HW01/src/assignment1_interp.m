figure; 
h = drawpolyline;
z = h.Position;
n = size(z, 1);
x = z(:,1);
y = z(:,2);

hold on;

fInterpPoly = @(x, y, xs) polyInterp(x, y, xs);

d0 = 0.1;
logd0 = log10(d0);    
lambda0 = 1e-12;

fInterpRBF = @(x, y, xs, d) rbfInterpConstTail(x, y, xs, d, lambda0);

function updateRBF(hrbf, x, y, xs, d, txt)
    set(hrbf, 'ydata', rbfInterpConstTail(x, y, xs, d, 1e-12));
    txt.String = sprintf('d = %.3e', d);
end

xs = (0:0.01:1)';
hpoly = plot(xs, fInterpPoly(x, y, xs), 'r', 'linewidth', 2);
hrbf = plot(xs, fInterpRBF(x, y, xs, d0), 'g', 'linewidth', 2);

% add slider for log10(d)
sld = uicontrol('Style', 'slider', ...
    'Min', -3, 'Max', 1, 'Value', logd0, ...   % d from 1e-3 to 10^1
    'Units','normalized', ...
    'Position', [0.2 0.01 0.6 0.05]); % bottom of figure
txt = uicontrol('Style','text', ...
    'Units','normalized', ...
    'Position',[0.82 0.01 0.15 0.05], ...
    'String', sprintf('d = %.3e', d0));

% update polynomial when points move
addlistener(h, 'MovingROI', @(src,evt) ...
    set(hpoly, 'ydata', fInterpPoly( ...
    evt.CurrentPosition(:,1), evt.CurrentPosition(:,2), xs)));

addlistener(h, 'MovingROI', @(src,evt) ...
    set(hrbf, 'ydata', fInterpRBF( ...
    evt.CurrentPosition(:,1), evt.CurrentPosition(:,2), xs, ...
    10^(sld.Value))));

addlistener(sld, 'ContinuousValueChange', @(src,evt) ...
    updateRBF(hrbf, h.Position(:,1), h.Position(:,2), xs, ...
    10^(src.Value), txt));

legend('Polynomial', 'RBF');
axis manual;

%%
function result = polyInterp(x, y, xs)
    % INPUT:
    %   x, y : Input points
    %   xs   : Independent variable needs to be interpolated
    % OUTPUT:
    %   result : Interpolation result

    n = length(x);

    % Build divided difference table
    divdiff = y(:);
    for j = 2:n
        divdiff(j:n) = ( ...
            divdiff(j:n) - divdiff(j-1:n-1)) ./ (x(j:n) - x(1:n-j+1));
    end

    % Evaluate interpolation at query points xs
    m = length(xs);
    result = zeros(m,1);
    for k = 1:m
        val = divdiff(n);
        for j = n-1:-1:1
            val = val * (xs(k)-x(j)) + divdiff(j);
        end
        result(k) = val;
    end
end

%%
function yhat = rbfInterpConstTail(x, y, xs, d, lambda)
    % Radial Basis Function interpolation with constant tail
    % and constraint sum(b_j) = 0.
    %
    % INPUT:
    %   x      : n x dim matrix of data point coordinates
    %   y      : n x 1 vector of target values
    %   xs     : m x dim matrix of query points
    %   d      : shape parameter for kernel (default 1e-1)
    %   lambda : small regularization parameter (default 1e-12)
    %
    % OUTPUT:
    %   yhat   : m x 1 vector of interpolated values
    %
    % Kernel: g(r) = 1 / (r^2 + d)
    
    if nargin < 5
        lambda = 1e-12;   % default regularization
    end

    if nargin < 4
        d = 0.1;         % default shape parameter
    end
    
    n = size(x, 1);
    m = size(xs,1);
    
    % Build kernel matrix A
    A = zeros(n,n);
    for i = 1:n
        for j = 1:n
            r2 = sum((x(i,:) - x(j,:)).^2);
            A(i,j) = 1 / (r2 + d);
        end
    end
    
    % Augmented system
    % [ A+lambda*I   1 ]
    % [   1^T        0 ]
    M = [A + lambda*eye(n), ones(n,1);
         ones(1,n),         0       ];
    
    rhs = [y; 0];
    
    % Solve system
    sol = M \ rhs;
    b = sol(1:n);
    c = sol(end);
    
    % Evaluate interpolation at query points xs
    yhat = zeros(m,1);
    for k = 1:m
        r2 = sum((xs(k,:) - x).^2, 2);
        G = 1 ./ (r2 + d);
        yhat(k) = G' * b + c;
    end
end
