function [fobj, fresidue, fsparsity] = getObjective3(A, S, X, beta)

E = A*S - X;
fresidue  = 0.5*sum(sum(E.^2));
fsparsity = beta*sum(sum(abs(S)));

fobj = fresidue + fsparsity;
