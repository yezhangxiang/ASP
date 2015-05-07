function [IX, solution] = balanceTraffic(IX, E, solution, gridSize, Tuenum, TrafficMap)

[row, col] = size(IX);
[M, ~] = size(solution);
blockTraffic = zeros(M, 1);
for i = 1 : M
    blockTraffic(i) = sum(sum(TrafficMap(IX == i)));
end
illBlock = 1: M;
illBlock = illBlock(blockTraffic>Tuenum);
if isempty(illBlock)
    return
end
for i = 1 : length(illBlock)
    currentIllBlock = illBlock(i);
    E1 = E(:, 1);
    E2 = E(:, 2);
    neighbor1 = E2(E1 == currentIllBlock);
    neighbor = [neighbor1; E1(E2 == currentIllBlock)];
    filter = ~ismember(neighbor, illBlock);
    neighbor = neighbor(filter);
    if isempty(neighbor)
        continue
    end
    bounds = struct('neighbor', cell(length(neighbor), 1), 'bound', [], 'traffic', []);
    totleSize = 0;
    totleTaffic = 0;
    for n_i = 1 : length(neighbor)
        bound = false(row, col);
        num = 0;
        traffic = 0;
        for r_i = 2 : row-1
            for c_i = 2 : col-1
                if IX(r_i, c_i) == currentIllBlock && ...
                    (IX(r_i+1, c_i+1) == neighbor(n_i) || ...
                     IX(r_i+1, c_i-1) == neighbor(n_i) || ...
                     IX(r_i-1, c_i+1) == neighbor(n_i) || ...
                     IX(r_i-1, c_i-1) == neighbor(n_i) )
                    bound(r_i, c_i) = true;
                    num = num+1;
                    traffic = traffic + TrafficMap(r_i, c_i);
                end
            end
        end
        bounds(n_i).neighbor = neighbor(n_i);
        bounds(n_i).bound = bound;
        bounds(n_i).num = num;
        bounds(n_i).traffic = traffic;
        totleSize = totleSize+num;
        totleTaffic = totleTaffic + traffic;
    end
    deltaV = [0, 0];
    averDisBetweenSite = 0;
    %if totleSize <= blockTraffic(currentIllBlock)-Tuenum
        for n_i = 1 : length(neighbor)
            IX(bounds(n_i).bound) = bounds(n_i).neighbor;
            deltaV_i = solution(currentIllBlock, 1:2)-solution(bounds(n_i).neighbor, 1:2);
            disBetweenSite = norm(deltaV_i);
            deltaV_i = deltaV_i./disBetweenSite.*0.5*gridSize;
            solution(bounds(n_i).neighbor, 1:2) = solution(bounds(n_i).neighbor, 1:2) + deltaV_i;
            deltaV = deltaV + deltaV_i;
            averDisBetweenSite = averDisBetweenSite + disBetweenSite;
        end
        deltaP = 0.5;
        if norm(deltaV) > gridSize/2
            deltaV = deltaV./norm(deltaV).*gridSize/2;

        end
        solution(currentIllBlock, 1:2) = solution(currentIllBlock, 1:2) + deltaV;
        solution(currentIllBlock, 3) = solution(currentIllBlock, 3)-(1-norm(deltaV)/gridSize*deltaP);
        averDisBetweenSite = averDisBetweenSite/length(neighbor);
    %end
    %solution(currentIllBlock, :) = solution(currentIllBlock, :) + deltaV;
end





