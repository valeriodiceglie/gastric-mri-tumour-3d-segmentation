function analyzeSeedFallback()
%ANALYZESEEDFALLBACK Effetto di seedFallback: il metodo finale contro la stessa configurazione senza.
%   Senza seedFallback, un seme che cade su un pixel sotto soglia non fa crescere nulla: 1 paziente
%   dev e 16 pazienti test risultavano vuoti nel primo passaggio sul test, fatto senza correzione.
%   Uscite in results/analyzeSeedFallback/: per_patient.csv, summary.csv
outDir = "results/analyzeSeedFallback"; [~, ~] = mkdir(outDir); tic;
[P1, Q] = methodConfig(); P0 = P1; P0.seedFallback = false;
rows = table(); summ = table();
for set = ["dev", "test"]
    T = readCases("manifest/" + set + "_cases.csv");
    T = sortrows(T(T.modality == "V", :), 'patient_id'); n = height(T);
    D0 = nan(n, 1); D1 = D0; E0 = false(n, 1); E1 = E0;
    for i = 1:n
        [I, S, sp] = readCase(T.img(i), T.seg(i));
        M0 = segmentCase(int16(I), S, sp, P0, Q); M1 = segmentCase(int16(I), S, sp, P1, Q);
        D0(i) = diceCoef(M0, S); D1(i) = diceCoef(M1, S); E0(i) = ~any(M0(:)); E1(i) = ~any(M1(:));
        if mod(i, 50) == 0, fprintf('%s %d/%d (%.0f s)\n', set, i, n, toc); end
    end
    rows = [rows; table(repmat(set, n, 1), T.patient_id, D0, D1, E0, E1, 'VariableNames', ...
        {'set', 'patient_id', 'dice3d_senza', 'dice3d_con', 'vuoto_senza', 'vuoto_con'})]; %#ok<AGROW>
    c0 = bootMedianCI(D0); c1 = bootMedianCI(D1); ch = abs(D1 - D0) > 1e-12;
    summ = [summ; table(set, n, median(D0), c0(1), c0(2), median(D1), c1(1), c1(2), mean(D0), mean(D1), ...
        100 * mean(D0 < 0.1), 100 * mean(D1 < 0.1), nnz(E0), nnz(E1), nnz(ch), wilcoxonSignedRank(D1 - D0), ...
        'VariableNames', {'set', 'n', 'mediana_senza', 'ci95_basso_senza', 'ci95_alto_senza', ...
        'mediana_con', 'ci95_basso_con', 'ci95_alto_con', 'media_senza', 'media_con', ...
        'fallimenti_pct_senza', 'fallimenti_pct_con', 'vuoti_senza', 'vuoti_con', 'pazienti_cambiati', 'p_wilcoxon'})]; %#ok<AGROW>
    fprintf('\n== %s (%d): senza %.3f (%.3f-%.3f) | con seedFallback %.3f (%.3f-%.3f) | vuoti %d -> %d | cambiati %d\n', ...
        set, n, median(D0), c0, median(D1), c1, nnz(E0), nnz(E1), nnz(ch));
    if any(ch), fprintf('   %s\n', strjoin(compose('%s %.3f->%.3f', T.patient_id(ch), D0(ch), D1(ch)), ', ')); end
end
disp(summ);
writetable(rows, fullfile(outDir, 'per_patient.csv'));
writetable(summ, fullfile(outDir, 'summary.csv'));
end
