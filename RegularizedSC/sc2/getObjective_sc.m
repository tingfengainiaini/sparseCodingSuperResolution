function [fobj, fresidue, freg] = getObjective_sc(A, S, X, beta)

E = double(A)*S - double(X);
fresidue = sum(sum(E.^2));
freg = beta*sum(sum(abs(S)));

fobj= fresidue + freg;


