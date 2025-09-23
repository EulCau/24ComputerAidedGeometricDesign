figure; 
h = impoly('closed', false);
z = h.getPosition;
n = size(z, 1);
x = z(:,1);
y = z(:,2);

hold on;

fInterpPoly = @(x, y, xs) polyInterp(x, y, xs);
fInterpRBF = @(x, y, xs) xs;

xs = (0:0.01:1)';
hpoly = plot(xs, fInterpPoly(x, y, xs), 'r', 'linewidth', 2);
hrbf = plot(xs, fInterpRBF(x, y, xs), 'g', 'linewidth', 2);

h.addNewPositionCallback(@(y) set(hpoly, 'ydata', fInterpPoly(y(:,1),y(:,2), xs)) );
h.addNewPositionCallback(@(y) set(hrbf, 'ydata', fInterpRBF(y(:,1),y(:,2), xs)) );

legend('Polynomial', 'RBF'); %axis manual;

%%

function yhat = polyInterp(x, y, xs)
    % Input:
    %   x, y : Input points
    %   xs   : Independent variable needs to be interpolated
    % output:
    %   yhat : Interpolation result

    n = length(x);
    % Make Vandermonde matrix
    V = zeros(n, n);
    for i = 1:n
        V(i, :) = x(i).^(0:n-1);
    end

    % Solve result
    alpha = V \ y;

    % Compute interpolation result
    Xs = zeros(length(xs), n);
    for i = 1:length(xs)
        Xs(i, :) = xs(i).^(0:n-1);
    end
    yhat = Xs * alpha;
end


%%
function [result] = rbf(n,x,y)


end

