function [objDF, objWav] = wObjective(b, x, lam, param)


objDF = param.A(x) - b;
objDF = objDF(:)'*objDF(:);
 
coef = param.U(x);
% objWav = sum(cell2mat(cellfun(@(x,y) (sum(sum(abs(x.*y)))), coef, lam, 'un',0)));
objWav = abs(coef(:).*lam(:));
objWav = sum(objWav(:));
