[x, t] = readObj('cathead.obj');

figure; subplot(121); trimesh(t, x(:,1), x(:,2), x(:,3), 'edgecolor', 'k'); axis equal; axis off; title('input mesh');

y = tutte_embed(x, t);

subplot(122); trimesh(t, y(:,1), y(:,2), y(:,1)*0, 'edgecolor', 'k'); axis equal; axis off; title('Tutte embedding'); view(2);


%% Tutte embedding algorithm
function y = tutte_embed(x, t)
    % 1. Identify Boundary Vertices
    % B is a list of boundary vertex IDs in order.
    B = findBoundary(x, t);
    nB = length(B);
    num_v = size(x, 1);
    
    % 2. Map Boundary to a Convex Shape (Unit Circle)
    % Using arc-length or simple angle distribution to fix boundary points.
    theta = linspace(0, 2*pi, nB+1)'; 
    theta(end) = []; % Remove duplicate end point
    
    % Boundary coordinates in 2D
    b_coords = [cos(theta), sin(theta)];
    
    % 3. Compute Weights (Floater / Shape Preserving)
    % We iterate over edges to build the weight matrix.
    % Edges in the triangulation:
    E = [t(:,1) t(:,2); t(:,2) t(:,3); t(:,3) t(:,1)];
    
    % Make edges unique (undirected) to avoid double counting during weight calc
    E = sort(E, 2);
    E = unique(E, 'rows');
    
    % Calculate Euclidean distance for every edge
    diffs = x(E(:,1), :) - x(E(:,2), :);
    dists = sqrt(sum(diffs.^2, 2));
    
    % SHAPE PRESERVING WEIGHTS (Inverse Distance)
    % w_ij = 1 / ||v_i - v_j||
    % (For standard Tutte/Uniform, you would use w = ones(size(dists)))
    weights = 1 ./ dists; 
    
    % 4. Construct the Laplacian/System Matrix
    % Build a symmetric sparse matrix W
    % We place weights at (i,j) and (j,i)
    W = sparse([E(:,1); E(:,2)], [E(:,2); E(:,1)], [weights; weights], num_v, num_v);
    
    % Compute the diagonal degree matrix D (row sums of W)
    D = spdiags(sum(W, 2), 0, num_v, num_v);
    
    % The Laplacian L = D - W
    L = D - W;
    
    % 5. Apply Boundary Conditions
    % For boundary vertices, we strictly enforce y_i = fixed_position.
    % We modify the L matrix and the RHS to solve: L * y = rhs
    
    rhs = zeros(num_v, 2);
    
    % Set boundary rows in L to Identity
    L(B, :) = 0; 
    
    % We use a linear index to set the diagonal elements for boundary rows to 1
    L(sub2ind(size(L), B, B)) = 1; 
    
    % Set the Right Hand Side (RHS) for boundary vertices
    rhs(B, :) = b_coords;
    
    % 6. Solve the Linear System
    % Solve L * y = rhs
    y = L \ rhs;
    
end
