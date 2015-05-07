function asp (fileName)

% global M RSCP x y row col A B Hue Trsrp Tsinr Tuenum TrafficMap w di stance
global distance
% fileName = 'input010.txt';
[A, B, Trsrp, Tsinr, Tuenum, w, M, TrafficMap] = readInput(fileName);
[row, col] = size(TrafficMap);
RSCP = zeros(M, row, col);
gridSize = 10;
x = gridSize/2 : gridSize : col*gridSize;
y = row*gridSize-gridSize/2 : -gridSize : 0;
Nue = sum(sum(TrafficMap));

x2 = gridSize/2 : gridSize : 2*col*gridSize;
y2 = 2*row*gridSize-gridSize/2 : -gridSize : 0;
distance = bsxfun(@plus, x2, (y2*1i).');
distance = abs(distance - (gridSize*col+gridSize*row*1i));
distance = A+B*log10(distance);

figure
imagesc(x, y, TrafficMap)
hold on
title('TrafficMap')
set(gca, 'YDir', 'normal');
colorbar
axis([min(x) max(x) min(y) max(y)])
axis equal

tStart = tic;

initPower = 15 .* ones(M,1);
%% 用kmeans算法给初值
[meshX, meshY]=meshgrid(x, y);
xd = reshape(meshX, col*row, 1);
yd = reshape(meshY, col*row, 1);
if M > 10
    kmeansTimes = 3;
else
    kmeansTimes = (10-M)*2-3;
end
opts = statset('MaxIter',200);
[IX, ctrs] = kmeans([xd yd], M, 'Replicates' ,kmeansTimes, 'Options' ,opts);
IX = reshape(IX, row, col);
kmeans_solution = [ctrs initPower];
kmeans_solution = regulariseSolution(kmeans_solution, x, y, 9:15, M);
[kmeans_fitness, RSCP] = gaFitnessfcn(kmeans_solution, 1:M, x, y, RSCP, TrafficMap, Nue, Trsrp, Tsinr, Tuenum, w, row, col);
figure('Name', fileName)
subplot(2, 3, 1)
imagesc(x, y, IX)
hold on
plot(kmeans_solution(: ,1) ,kmeans_solution(: ,2) ,'kx', 'MarkerSize' ,12, 'LineWidth', 2)
title(['Phase 1: '  num2str(kmeans_fitness)])
set(gca, 'YDir',' normal' ) ;
axis equal
disp(['after k-means fitness is ' num2str(kmeans_fitness)])


DT = delaunayTriangulation(ctrs);
E = edges(DT);
solution = kmeans_solution;
for ii = 1 : 30
    [IX, solution] = balanceTraffic(IX, E, solution, gridSize, Tuenum, TrafficMap);
end

traffic_solution = solution;
traffic_solution = regulariseSolution(traffic_solution, x, y, 9:15, M);

traffic_solution2 = solution;
tmpX = repmat(x, row, 1);
tmpY = repmat(y', 1, col);
for i = 1 : M
    tmpIndex = IX == i;
    traffic_solution2(i, 1) = mean(tmpX(tmpIndex));
    traffic_solution2(i, 2) = mean(tmpY(tmpIndex));
end
traffic_solution2 = regulariseSolution(traffic_solution2, x, y, 9:15, M);

[fitness1, RSCP] = gaFitnessfcn(traffic_solution, 1:M, x, y, RSCP, TrafficMap, Nue, Trsrp, Tsinr, Tuenum, w, row, col);
[fitness2, RSCP2] = gaFitnessfcn(traffic_solution2, 1:M, x, y, RSCP, TrafficMap, Nue, Trsrp, Tsinr, Tuenum, w, row, col);
if fitness2 > fitness1
    traffic_solution = traffic_solution2;
    RSCP = RSCP2;
end

% 显示traffic_solution
[besTrsrp, bestCI] = max(RSCP);
besTrsrp = reshape(besTrsrp, row, col);
bestCI = reshape(bestCI, row, col);
SINR = calSINR(besTrsrp, RSCP);
traffic_fitness = calFitness(w, M, besTrsrp, bestCI, SINR, Trsrp, Tsinr, Tuenum, TrafficMap, Nue, row, col);

disp(['after traffic average fitness is ' num2str(traffic_fitness)])

subplot(2, 3, 2)
imagesc(x, y, IX)
hold on
plot(traffic_solution(:, 1), traffic_solution(:, 2), 'kx', 'MarkerSize', 12, 'LineWidth', 2)
title('Phase 2: prediction')
set(gca, 'YDir', 'normal');
axis equal

subplot(2, 3, 3)
imagesc(x, y, bestCI)
hold on
plot(traffic_solution(:, 1), traffic_solution(:, 2), 'kx', 'MarkerSize', 12, 'LineWidth', 2)
title(['Phase 2: '  num2str(traffic_fitness)])
set(gca, 'YDir', 'normal');
axis equal

subplot(2, 3, 4)
Fgrid = besTrsrp >= Trsrp & SINR >= Tsinr;
imagesc(x, y, Fgrid)
hold on
plot(traffic_solution(:, 1), traffic_solution(:, 2), 'kx', 'MarkerSize', 12, 'LineWidth', 2)
title('Phase 2: problem')
set(gca, 'YDir', 'normal');
axis equal

%%寻优
xyStep = gridSize;
xySearchScope = round(sqrt(row*col/(M*pi))/2);
pSearchScope = 5;
groupSize = 3; %2小区一组
groupNum = ceil(M/groupSize);
if traffic_fitness > kmeans_fitness
    solution = traffic_solution;
else
    solution = kmeans_solution;
    for i = 1 : M
        RSCP(i, :, :) = sigleCalRSCP(x, y, solution(i, 1), solution(i, 2), solution(i, 3));
    end
end

%% 分组
[besTrsrp, bestCI] = max(RSCP);
besTrsrp = reshape(besTrsrp, row, col);
bestCI = reshape(bestCI, row, col);
SINR = calSINR(besTrsrp, RSCP);
tmpRSCP = reshape(RSCP, M, row*col);
groupIX = kmeans(tmpRSCP, groupNum);
% 组排序
group_score = zeros(groupNum, 1);
problemGrid = besTrsrp<Trsrp | SINR <Tsinr;
for group_i = 1 : groupNum
    cell_i = groupIX == group_i;
    tmpCell = 1 :M;
    cell_i = tmpCell(cell_i);

    group_score(group_i) = sum(sum(ismember(bestCI, cell_i) & problemGrid));
end
[~, groupSort] = sort(group_score, 'descend');

DT = delaunayTriangulation(solution(:, 1:2));
E = edges(DT);
E1 = E(:, 1);
E2 = E(:, 2);
calTime = 0;
isTimeOut = false;
tOptStart = tic;

for mainIter = 1 : max(floor(log2(xySearchScope)), 5)


for group_i = 1 : groupNum
    % cell_i = TRI(group_i, :);
    cell_i = groupIX == groupSort(group_i);
    currentGroupSize = sum(cell_i);
    cell_x = solution(cell_i, 1);
    cell_y = solution(cell_i, 2);
    cell_p = solution(cell_i, 3);
    x_lower = cell_x - xySearchScope * gridSize;
    x_lower(x_lower < 5) = 5;
    x_upper = cell_x + xySearchScope * gridSize;
    x_upper(x_upper>col*gridSize-5) = col*gridSize-5;
    y_lower = cell_y - xySearchScope * gridSize;
    y_lower(y_lower < 5) = 5;
    y_upper = cell_y + xySearchScope * gridSize;
    y_upper(y_upper>row*gridSize-5) =row*gridSize-5;
    p_lower = cell_p - pSearchScope;
    p_lower(p_lower<9)=9;
    p_upper = cell_p + pSearchScope;
    p_upper(p_upper>15)=15;

    UB = [x_upper; y_upper; p_upper];
    LB = [x_lower; y_lower; p_lower];
    step = [xyStep.*ones(2*currentGroupSize, 1); ones(currentGroupSize, 1)];
    tmpCell = 1 : M;
    cell_i = tmpCell(cell_i);
    % 构造补丁
    tmpD = [];
    for tmpNeiI = 1 : length(cell_i)
        current_cell_i = cell_i(tmpNeiI);
        tmpNei1 = E2(E1 == current_cell_i);
        tmpNei2 = E1(E2 == current_cell_i);
        tmpNei = [tmpNei1; tmpNei2];
        tmpNeiX = solution(tmpNei, 1);
        tmpNeiY = solution(tmpNei, 2);
        D = abs(tmpNeiX + tmpNeiY.*1i - (solution(current_cell_i, 1)+solution(current_cell_i, 2).*1i));
        tmpD = [tmpD; D];
    end
    effectRadio = max(tmpD)*0.75;
    effectAreaXU = min([gridSize*col, max(solution(cell_i, 1))+effectRadio]);
    effectAreaXD = max([0, min(solution(cell_i, 1))-effectRadio]);
    effectAreaYU = min([gridSize*row, max(solution(cell_i, 2))+effectRadio]);
    effectAreaYD = max([0, min(solution(cell_i, 2))-effectRadio]);
    localX = round(effectAreaXD/gridSize)*gridSize + gridSize/2 : gridSize : effectAreaXU;
    localY = round(effectAreaYU/gridSize)*gridSize - gridSize/2 : -gridSize : effectAreaYD;
    localRSCP = RSCP(:, ismember(y, localY), ismember(x, localX));
    localTraffic = TrafficMap(ismember(y, localY), ismember(x, localX));





    dimCI = [cell_i';cell_i';cell_i'];
    dimPara = struct('scope', cell(3*currentGroupSize, 1), 'cellId', []);
    for i = 1 : 3*currentGroupSize
        dimPara(i).scope = LB(i) : step(i) : UB(i);
        dimPara(i).cellId = dimCI(i);
    end

    [maxScore, swarm_best, claNum] = optimiseGPSO(solution, dimPara, localX, localY, localRSCP, localTraffic, Nue, Trsrp, Tsinr, Tuenum, w, row, col);
    solution = para2solution(dimPara, solution, swarm_best);
    calTime = calTime+1;
    % 判断是否提前退出
    tElapsed = toc(tStart);
    tOptElapsed = toc(tOptStart);
    if (300-tElapsed)<tOptElapsed/calTime*3
        isTimeOut = true;
        break
    end

end
if isTimeOut
    break
end
% 每次大迭代缩小搜索范围
if xySearchScope > 2
    xySearchScope = ceil(xySearchScope/2);
end
if pSearchScope > 2
    pSearchScope = ceil(pSearchScope/2);
end
end
best_solution = solution;
for i = 1 : M
    RSCP(i, :, :) = sigleCalRSCP(x, y, solution(i, 1), solution(i, 2), solution(i, 3));
end
[besTrsrp, bestCI] = max(RSCP);
besTrsrp = reshape(besTrsrp, row, col);
bestCI = reshape(bestCI, row, col);
SINR = calSINR(besTrsrp, RSCP);
fitness = calFitness(w, M, besTrsrp, bestCI, SINR, Trsrp, Tsinr, Tuenum, TrafficMap, Nue, row, col);

subplot(2, 3, 5)
imagesc(x, y, bestCI)
hold on
plot(best_solution(:, 1), best_solution(:, 2), 'kx', 'MarkerSize', 12, 'LineWidth', 2)
title(['Phase 3: ' num2str(fitness)]);
set(gca, 'YDir', 'normal');
axis equal

subplot(2, 3, 6)
Fgrid = besTrsrp >= Trsrp & SINR >= Tsinr;
imagesc(x, y, Fgrid)
hold on
plot(best_solution(:, 1), best_solution(:, 2), 'kx', 'MarkerSize', 12, 'LineWidth', 2)
title('Phase 3: problem');
set(gca, 'YDir', 'normal');
axis equal
disp(['best fitness is ' num2str(fitness)])

final_solution = [(1:M)' best_solution];

final_solution(:, 2) = final_solution(:, 2)-gridSize/2;
final_solution(:, 3) = final_solution(:, 3)+gridSize/2;
outputFile = ['out' fileName(3:length(fileName))];
dlmwrite(outputFile, fitness, 'newline', 'pc');
dlmwrite(outputFile, final_solution, '-append', 'delimiter', '\t', 'newline', 'pc');

tElapsed = toc(tStart);
disp(['deal ' fileName ' calculate ' num2str(calTime) ' times, take ' num2str(tElapsed) 'seconds']);
