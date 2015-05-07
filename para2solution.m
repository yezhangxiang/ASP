function [solution, changeCI] = para2solution(dimPara, solution, cellPara)
dimNum = length(dimPara);
groupSize = dimNum/3;
changeCI = [];
for i = 1 : length(dimPara)
    solution(dimPara(i).cellId, ceil(i/groupSize)) = cellPara(i);
    if ~ismember(dimPara(i).cellId, changeCI)
        changeCI = [changeCI; dimPara(i).cellId];
    end
end
