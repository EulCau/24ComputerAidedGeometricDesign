%% Assignment 6 - Problem 1
% ------------------------------------------------------------
% Goal:
% (1) Use de Casteljau algorithm to evaluate sample points at u={2.5,3,3.5}
% (2) Subdivide the cubic Bezier curve at u=3 and u=3.5, and plot
%     all corresponding control polygons and curve segments
% ------------------------------------------------------------

clear; clc; close all;

%% --- Bézier control points from Problem 1(1)
% These are derived analytically for interval [2,4]
P = [2 0; 0 2; 3 4; 4 1];  % [x y]
a = 2; b = 4;               % parameter interval

%% --- Helper: de Casteljau evaluation
function Q = deCasteljau(P, t)
    n = size(P,1);
    for r = 1:n-1
        for i = 1:n-r
            P(i,:) = (1-t)*P(i,:) + t*P(i+1,:);
        end
    end
    Q = P(1,:);
end

%% --- Helper: subdivision (returns left/right control polygons)
function [left,right] = subdivide(P,t)
    n = size(P,1);
    tmp = P;
    left = zeros(n,2); right = zeros(n,2);
    left(1,:) = tmp(1,:);
    right(n,:) = tmp(n,:);
    for r = 1:n-1
        for i = 1:n-r
            tmp(i,:) = (1-t)*tmp(i,:) + t*tmp(i+1,:);
        end
        left(r+1,:) = tmp(1,:);
        right(n-r,:) = tmp(n-r,:);
    end
end

%% --- (2) Evaluate sample points at u = {2.5, 3, 3.5}
u = [2.5, 3.0, 3.5];
t = (u - a) / (b - a);   % map u to [0,1]

Q = zeros(length(t),2);
for k = 1:length(t)
    Q(k,:) = deCasteljau(P, t(k));
end
disp('Sample points [x y] at u={2.5,3,3.5}:');
disp(Q);

%% --- (3) Subdivide curve at u=3 and u=3.5
t1 = (3 - a)/(b - a);    % 0.5
t2 = (3.5 - a)/(b - a);  % 0.75

% First subdivision at u=3
[L1, R1] = subdivide(P, t1);

% Subdivide right half at u=3.5 (local t=0.5 within [3,4])
[L2, R2] = subdivide(R1, 0.5);

%% --- Draw everything
figure; hold on; axis equal; grid on;
title('Problem 1: de Casteljau Sampling and Subdivision');
xlabel('x'); ylabel('y');

% --- Plot full original Bezier curve
tt = linspace(0,1,200);
curve = zeros(length(tt),2);
for k = 1:length(tt)
    curve(k,:) = deCasteljau(P, tt(k));
end
plot(P(:,1),P(:,2),'k--o','LineWidth',1.2,'DisplayName','Original control polygon');
plot(curve(:,1),curve(:,2),'k','LineWidth',1.5,'DisplayName','Original Bezier curve');

% --- Plot subdivision polygons and curves
colors = {'r','g','b'};
segments = {L1,L2,R2};
labels = {'Left segment [2,3]','Middle [3,3.5]','Right [3.5,4]'};

for k = 1:3
    % control polygon
    plot(segments{k}(:,1),segments{k}(:,2),[colors{k} 'o--'],'LineWidth',1.3,'DisplayName',[labels{k} ' polygon']);
    % curve
    tt = linspace(0,1,100);
    subCurve = zeros(length(tt),2);
    for i = 1:length(tt)
        subCurve(i,:) = deCasteljau(segments{k}, tt(i));
    end
    plot(subCurve(:,1),subCurve(:,2),colors{k},'LineWidth',1.8,'DisplayName',[labels{k} ' curve']);
end

% --- Plot sample points from part (2)
plot(Q(:,1),Q(:,2),'ms','MarkerFaceColor','m','DisplayName','Sample points u={2.5,3,3.5}');

legend('Location','northwest');
