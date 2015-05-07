function fitness = calFitness(w, M, bestRSCP, bestCI, SINR, Trscp, Tsinr, Tuenum, TrafficMap, Nue, row, col)
Fgrid = bestRSCP >= Trscp & SINR >= Tsinr;
fitnessArea = sum(sum(Fgrid))/(row*col);
FBNUM = zeros(1, M);
for i = 1 : M
    Fue = (bestCI == i & Fgrid);
    FBNUM(i) = min(sum(sum(TrafficMap(Fue))), Tuenum);
end

fitnessUE = sum(FBNUM)/Nue;
fitness = w*fitnessArea + (1-w)*fitnessUE;
