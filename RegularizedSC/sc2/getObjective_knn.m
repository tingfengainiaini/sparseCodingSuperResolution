function [fobj, fresidue, freg] = getObjective_knn(A, S, X, beta)

E = double(A)*S - double(X);

fresidue = sum(sum(E.^2));
freg = beta*sum(sum(S.^2));

fobj= fresidue + freg;


