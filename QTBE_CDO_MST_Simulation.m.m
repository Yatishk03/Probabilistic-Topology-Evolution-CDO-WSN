% -------------------------------------------------------------------------
% QTBE + CDO + Kruskal MST Simulation in WSN (Multiple Rounds)
% Metrics: Node Death Rate, Network Lifetime, Energy Consumption per Round
% -------------------------------------------------------------------------

clc; clear; close all;

% Parameters
numNodes = 50;
fieldSize = 100;
numCHs = 6;
BS = [120, 50];  % Base Station out of field
R = 25;          % Tunneling range parameter
initialEnergy = 1;
energyThreshold = 0.01;

% Deployment
nodes = rand(numNodes, 2) * fieldSize;
energy = ones(1, numNodes) * initialEnergy;
colors = lines(numCHs);  % For cluster coloring

% Stats initialization
deathHistory = [];
energyHistory = [];
round = 1;

while any(energy > energyThreshold)
    % QTBE-based Probabilistic Topology
    distMatrix = squareform(pdist(nodes));
    T = exp(-distMatrix / R);
    G = zeros(numNodes);
    for i = 1:numNodes
        for j = i+1:numNodes
            if rand < T(i,j)
                G(i,j) = 1; G(j,i) = 1;
            end
        end
    end

    % CDO - Cluster Head Selection
    connectivity = sum(G, 2)';
    avgDist = mean(distMatrix, 2)';
    alpha = 0.4; beta = 0.3; gamma = 0.3;
    utility = alpha * energy + beta * connectivity - gamma * avgDist;
    utility(energy <= energyThreshold) = -inf;
    [~, sortedIdx] = sort(utility, 'descend');
    CHs = sortedIdx(1:numCHs);

    % Cluster Formation
    clusters = zeros(1, numNodes);
    for i = 1:numNodes
        if energy(i) <= energyThreshold, continue; end
        dists = vecnorm(nodes(CHs,:) - nodes(i,:), 2, 2);
        [~, idx] = min(dists);
        clusters(i) = CHs(idx);
    end

    % Kruskal MST on CHs + BS
    CH_positions = nodes(CHs, :);
    CH_plus_BS = [CH_positions; BS];
    n = size(CH_plus_BS, 1);
    edgeList = [];
    for i = 1:n
        for j = i+1:n
            d = norm(CH_plus_BS(i,:) - CH_plus_BS(j,:));
            edgeList = [edgeList; i, j, d];
        end
    end

    edgeList = sortrows(edgeList, 3);
    parent = 1:n;
    mstEdges = [];
    for i = 1:size(edgeList,1)
        u = edgeList(i,1); v = edgeList(i,2);
        pu = find_root(u, parent); pv = find_root(v, parent);
        if pu ~= pv
            parent(pv) = pu;
            mstEdges = [mstEdges; u, v];
        end
    end

    % Visualization
    figure(1); clf; hold on; axis equal;
    title(['Round ', num2str(round)]);
    [edge_i, edge_j] = find(G);
    for k = 1:length(edge_i)
        plot([nodes(edge_i(k),1), nodes(edge_j(k),1)], ...
             [nodes(edge_i(k),2), nodes(edge_j(k),2)], 'Color', [0.85 0.85 0.85]);
    end

    for i = 1:numNodes
        if energy(i) <= energyThreshold, continue; end
        if any(CHs == i)
            scatter(nodes(i,1), nodes(i,2), 100, 'r', 'filled');
        else
            chIndex = find(CHs == clusters(i));
            scatter(nodes(i,1), nodes(i,2), 50, colors(chIndex,:), 'filled');
        end
    end

    scatter(BS(1), BS(2), 150, 'g', 'filled');
    for i = 1:size(mstEdges,1)
        a = mstEdges(i,1); b = mstEdges(i,2);
        plot([CH_plus_BS(a,1), CH_plus_BS(b,1)], ...
             [CH_plus_BS(a,2), CH_plus_BS(b,2)], 'k-', 'LineWidth', 2);
    end

    legend('QTBE Links','Sensor Nodes','Cluster Heads','Base Station','MST Links');
    xlabel('X'); ylabel('Y'); grid on;
    pause(0.1);

    % Energy Consumption
    E_tx = 0.001;
    totalEnergy = 0;
    for i = 1:numNodes
        if energy(i) <= energyThreshold || clusters(i) == 0, continue; end
        d = norm(nodes(i,:) - nodes(clusters(i),:));
        consumption = E_tx * d;
        energy(i) = energy(i) - consumption;
        totalEnergy = totalEnergy + consumption;
    end
    for ch = CHs
        if energy(ch) <= energyThreshold, continue; end
        d = norm(nodes(ch,:) - BS);
        consumption = E_tx * d;
        energy(ch) = energy(ch) - consumption;
        totalEnergy = totalEnergy + consumption;
    end

    % Update Stats
    deadNodes = sum(energy <= energyThreshold);
    deathHistory = [deathHistory, deadNodes];
    energyHistory = [energyHistory, totalEnergy];
    round = round + 1;
end

% Final Statistics
figure;
subplot(2,1,1);
plot(1:length(deathHistory), deathHistory, '-r', 'LineWidth', 2);
title('Node Death Over Time');
xlabel('Round'); ylabel('Number of Dead Nodes'); grid on;

subplot(2,1,2);
plot(1:length(energyHistory), energyHistory, '-b', 'LineWidth', 2);
title('Energy Consumption Per Round');
xlabel('Round'); ylabel('Energy Used'); grid on;

fprintf('Network Lifetime: %d rounds\n', round - 1);
