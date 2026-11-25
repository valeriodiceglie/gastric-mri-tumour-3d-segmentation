function d = diceCoef(A, B)
%DICECOEF Coefficiente di Dice tra due maschere logiche: 2|A∩B| / (|A|+|B|).
    A = logical(A); B = logical(B);
    s = nnz(A) + nnz(B);
    if s == 0, d = 1; else, d = 2*nnz(A & B)/s; end
end
