function [maxScore, swarm_best, calNum] = optimiseGPSO(solution, dimPara, x, y, RSCP, TrafficMap, Nue, Trsrp, Tsinr, Tuenum, w, row, col)










dimNum = length(dimPara);
groupSize = dimNum/3;
calNum = 0;
PsoPara.stop_criterion = 1;
PsoPara.iteration_num = 200;
PsoPara.c1 = 0.8;
PsoPara.c2 = 0.75;
PsoPara.w = 0.15;
PsoPara.swarm_size = 20;

swarm = zeros(PsoPara.swarm_size, dimNum);
orignalPartical = zeros(1, dimNum);
for i  = 1 : dimNum
    orignalPartical(i) = solution(dimPara(i).cellId, ceil(i/groupSize));
    tmpScope = dimPara(i).scope;
    scoreSize = length(tmpScope);
    tmpIndex = randi(scoreSize, PsoPara.swarm_size, 1);
    swarm(:, i) = tmpScope(tmpIndex);
end

swarm(1, :) = orignalPartical;
swarm_tested = zeros(PsoPara.swarm_size * PsoPara.iteration_num, dimNum);

particle_best_fitness = zeros(1, PsoPara.swarm_size);
swarm_best_fitness = 0;
for i = 1: PsoPara.swarm_size
    cellPara = swarm(i, :);
    [solution, changeCI] = para2solution(dimPara, solution, cellPara);
    score = gaFitnessfcn(solution, changeCI, x, y, RSCP, TrafficMap, Nue, Trsrp, Tsinr, Tuenum, w, row, col);
    calNum = calNum + 1;
    swarm_tested(calNum, :) = cellPara;
    particle_best_fitness(i) = score;
    if score > swarm_best_fitness
        swarm_best = cellPara;
        swarm_best_fitness = score;
    end
end

particle_best=swarm;

for iter = 1 : PsoPara.iteration_num

    % 种群更新
    % 1�?以概率c1*rand()与个体最优解进行交叉操作
    particleBestCrossIndex = rand(PsoPara.swarm_size, dimNum) < PsoPara.c1;
    swarm(particleBestCrossIndex) = particle_best(particleBestCrossIndex);
    % 2�?以概率c1*rand()与个全局�?��解进行交叉操�?
    swarmBestCrossIndex = rand(PsoPara.swarm_size, dimNum) < PsoPara.c2;
    tmpSwarmBest = repmat(swarm_best, PsoPara.swarm_size, 1);
    swarm(swarmBestCrossIndex) = tmpSwarmBest(swarmBestCrossIndex);
    % 3�?以概率w 进行变异操作
    variationSwarm = zeros(PsoPara.swarm_size, dimNum);
    for i = 1: dimNum
        tmpScope = dimPara(i).scope;
        scoreSize = length(tmpScope);
        tmpIndex = randi(scoreSize, PsoPara.swarm_size, 1);
        variationSwarm(:, i) = tmpScope(tmpIndex);
    end
    variationIndex = rand(PsoPara.swarm_size, dimNum) < PsoPara.w;
    swarm(variationIndex) = variationSwarm(variationIndex);

    for i = 1 : PsoPara.swarm_size
        cellPara = swarm(i, :);
        if ismember(cellPara, swarm_tested, 'rows')
            continue
        end
        [solution, changeCI] = para2solution(dimPara, solution, cellPara);
        score = gaFitnessfcn(solution, changeCI, x, y, RSCP, TrafficMap, Nue, Trsrp, Tsinr, Tuenum, w, row, col);

        calNum = calNum + 1;
        swarm_tested(calNum, :) = cellPara;

        % 个体�?��解更�?
        if score > particle_best_fitness(i)
            particle_best(i, :) = cellPara;
            particle_best_fitness(i) = score;
        end

        % 群体�?��解更�?
        if score > swarm_best_fitness
            swarm_best = cellPara;
            swarm_best_fitness = score;
        end
    end
    all_best_swarm_fitness(iter) = swarm_best_fitness;
    all_best_particle(iter, :) = swarm_best;
    % 如果连续5次迭代都没有增长就�?�?
    if (PsoPara.stop_criterion > 0 && iter > 5 && abs(all_best_swarm_fitness(iter)-all_best_swarm_fitness(iter-5)) <= 0)
        break;
    end
end
maxScore = swarm_best_fitness;
swarm_tested = swarm_tested(1 : calNum, :);





