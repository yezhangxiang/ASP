function [fitness, RSCP] = gaFitnessfcn(answer, changeCI, x, y, RSCP, TrafficMap, Nue, Trscp, Tsinr, Tuenum, w, row, col)

[M, m, n] = size(RSCP);
answer = reshape(answer, M, 3);
for i = 1 : length(changeCI)
    ci = changeCI(i);
    RSCP(ci, :, :) = sigleCalRSCP(x, y, answer(ci, 1), answer(ci, 2), answer(ci, 3));
end
[bestRSCP, bestCI] = max(RSCP);
bestRSCP = reshape(bestRSCP, m, n);
bestCI = reshape(bestCI, m, n);
SINR = calSINR(bestRSCP, RSCP);
fitness = calFitness(w, M, bestRSCP, bestCI, SINR, Trscp, Tsinr, Tuenum, TrafficMap, Nue, row, col);

