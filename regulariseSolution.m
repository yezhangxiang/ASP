function solution = regulariseSolution(solution, x, y, p, M)
% 规则化solution
currentDim = solution(:, 1);
currentScope = x;
tmpScope = repmat(currentScope, M, 1);
tmp = tmpScope - repmat(currentDim, 1, length(currentScope));
[~, I] = min(abs(tmp), [], 2);
solution(:, 1) = currentScope(I);

currentDim = solution(:, 2);
currentScope = y;
tmpScope = repmat(currentScope, M, 1);
tmp = tmpScope - repmat(currentDim, 1, length(currentScope));
[~, I] = min(abs(tmp), [], 2);
solution(:, 2) = currentScope(I);

currentDim = solution(:, 3);
currentScope = p;
tmpScope = repmat(currentScope, M, 1);
tmp = tmpScope - repmat(currentDim, 1, length(currentScope));
[~, I] = min(abs(tmp), [], 2);
solution(:, 3) = currentScope(I);
