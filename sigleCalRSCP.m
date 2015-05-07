function sigleRSCP = sigleCalRSCP(x, y, cellX, cellY, Pb)
global distance
[row, col] = size(distance);
row = row/2;
col = col/2;

m = length(y);
n = length(x);
gridSize = 10;
deltaX = round((cellX-x(1))/gridSize);
deltaY = round((y(1)-cellY)/gridSize);

distanceTmp = distance(row-deltaY+1:row-deltaY+m, col-deltaX+1:col-deltaX+n);
sigleRSCP = Pb - distanceTmp;
