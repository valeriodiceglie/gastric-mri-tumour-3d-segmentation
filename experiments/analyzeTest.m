function analyzeTest()
%ANALYZETEST Valutazione del metodo finale (segmentCase) sul set test, 474 pazienti.
%   Nessun parametro viene modificato. Calcola:
%     - Dice per volume con IC 95% della mediana, quota di fallimenti (Dice < 0,1)
%     - precisione e recall per voxel: quanto del risultato e' tumore, quanto tumore e' trovato
%     - Dice per slice secondo la posizione nel tumore (0% = prima slice, 100% = ultima)
%     - sottogruppi: volume del tumore, spessore delle slice, pixel, numero di semi, test originale
%     - volume stimato contro volume dell'esperto (errore relativo, correlazione di Spearman)
%   Uscite in results/analyzeTest/: per_patient.csv, dice_by_position.csv, subgroups.csv,
%   volumes.csv, test_distribution.png, volumes_precision_recall.png, cases_best_median_worst.png
T = readCases("manifest/test_cases.csv");
T = sortrows(T(T.modality == "V", :), 'patient_id'); n = height(T); ids = T.patient_id;
outDir = "results/analyzeTest"; [~, ~] = mkdir(outDir); tic;

%% Metodo finale su ogni paziente
D3 = nan(n, 1); SPX = D3; SPZ = D3; VE = D3; VP = D3; NSEED = D3; FP = D3; REC = D3; PREC = D3; SENS = D3;
SL = zeros(0, 3);                                       % [paziente, posizione nel tumore, Dice della slice]
for i = 1:n
    [I, S, sp] = readCase(T.img(i), T.seg(i));
    [M3, ~, sd] = segmentCase(int16(I), S, sp); NSEED(i) = size(sd, 1);
    vm = volMetrics(M3, S); D3(i) = vm(1); REC(i) = vm(2); FP(i) = vm(3);
    SPX(i) = sp(1); SPZ(i) = sp(3);
    VE(i) = nnz(S) * prod(sp) / 1000; VP(i) = nnz(M3) * prod(sp) / 1000;      % cm^3
    PREC(i) = nnz(M3 & S) / nnz(M3); SENS(i) = nnz(M3 & S) / nnz(S);         % precisione NaN se M3 e' vuota
    ks = find(squeeze(any(S, [1 2])));
    for k = ks'
        pos = 0.5; if numel(ks) > 1, pos = (k - ks(1)) / (ks(end) - ks(1)); end
        SL(end+1, :) = [i, pos, diceCoef(M3(:, :, k), S(:, :, k))]; %#ok<AGROW>
    end
    if mod(i, 25) == 0, fprintf('%d/%d pazienti (%.0f s)\n', i, n, toc); end
end

%% Risultato principale
ci = bootMedianCI(D3);
fprintf('\n== Test: %d pazienti ==\nDice per volume: mediana %.3f (IC95%% %.3f-%.3f), media %.3f, q25 %.3f, q75 %.3f\n', ...
    n, median(D3), ci(1), ci(2), mean(D3), prctile(D3, 25), prctile(D3, 75));
fprintf('Fallimenti (Dice < 0.1): %.1f%%  | slice del tumore trovate (mediana) %.2f | slice senza tumore segmentate (media) %.2f\n', ...
    100 * mean(D3 < 0.1), median(REC), mean(FP));
fprintf('Precisione per voxel: mediana %.3f (q25 %.3f, q75 %.3f; %d risultati vuoti esclusi) | recall per voxel: mediana %.3f (q25 %.3f, q75 %.3f)\n', ...
    median(PREC, 'omitnan'), prctile(PREC, 25), prctile(PREC, 75), nnz(isnan(PREC)), median(SENS), prctile(SENS, 25), prctile(SENS, 75));

%% Dice per slice secondo la posizione nel tumore
edges = [0 0.2 0.4 0.6 0.8 1.0001]; lab = ["0-20%"; "20-40%"; "40-60%"; "60-80%"; "80-100%"];
b = discretize(SL(:, 2), edges);
Spos = table(lab, accumarray(b, SL(:, 3), [5 1], @median), accumarray(b, 1, [5 1]), ...
    'VariableNames', {'posizione_nel_tumore', 'dice_slice_mediana', 'n_slice'});
fprintf('\n== Dice per slice secondo la posizione (0 = prima slice del tumore, 100 = ultima) ==\n'); disp(Spos);
writetable(Spos, fullfile(outDir, 'dice_by_position.csv'));

%% Sottogruppi
sf = readCases('manifest/split_final.csv');
t40 = ismember(ids, sf.patient_id(lower(sf.test40_originale) == "true" | sf.test40_originale == "1"));
tv = prctile(VE, [100/3 200/3]);
grp = {"Volume tumore: piccolo", VE <= tv(1); "Volume tumore: medio", VE > tv(1) & VE <= tv(2); ...
       "Volume tumore: grande", VE > tv(2); "Spessore slice <= 2.5 mm", SPZ <= 2.5 + 1e-6; ...
       "Spessore slice > 2.5 mm", SPZ > 2.5 + 1e-6; "Pixel < 0.8 mm", SPX < 0.8; "Pixel >= 0.8 mm", SPX >= 0.8; ...
       "Un solo seme", NSEED == 1; "Piu' semi", NSEED > 1; "Test originale (40)", t40; "Resto del test", ~t40};
Sg = table('Size', [size(grp, 1) 6], 'VariableTypes', {'string', 'double', 'double', 'double', 'double', 'double'}, ...
    'VariableNames', {'sottogruppo', 'n', 'dice3d_median', 'ci95_basso', 'ci95_alto', 'quota_fallimenti'});
for g = 1:size(grp, 1)
    x = D3(grp{g, 2}); c = bootMedianCI(x);
    Sg(g, :) = {grp{g, 1}, numel(x), median(x), c(1), c(2), mean(x < 0.1)};
end
fprintf('\n== Sottogruppi (terzili di volume: %.1f e %.1f cm3) ==\n', tv(1), tv(2)); disp(Sg);
writetable(Sg, fullfile(outDir, 'subgroups.csv'));
writetable(table(ids, D3, PREC, SENS, VE, SPX, SPZ, NSEED, REC, FP, t40, 'VariableNames', {'patient_id', 'dice3d', ...
    'precisione_voxel', 'recall_voxel', 'volume_cm3', 'pixel_mm', 'slice_mm', 'n_semi', 'slice_trovate', 'slice_fp', 'test40'}), fullfile(outDir, 'per_patient.csv'));

%% Volumi: stimato contro esperto
rel = (VP - VE) ./ VE;
rs = corrcoef(ranks(VE), ranks(VP)); rs = rs(1, 2);                   % Spearman = Pearson sui ranghi
fprintf('\n== Volumi ==\nVolume esperto mediano %.1f cm3, stimato %.1f cm3\n', median(VE), median(VP));
fprintf('Errore relativo con segno: mediana %+.0f%% (q25 %+.0f%%, q75 %+.0f%%)\n', 100 * median(rel), ...
    100 * prctile(rel, 25), 100 * prctile(rel, 75));
fprintf('Errore relativo assoluto: mediana %.0f%%;  correlazione di Spearman %.2f\n', 100 * median(abs(rel)), rs);
fprintf('Pazienti sottostimati: %.0f%%\n', 100 * mean(rel < 0));
writetable(table(ids, VE, VP, rel, D3, 'VariableNames', {'patient_id', 'volume_esperto_cm3', 'volume_stimato_cm3', ...
    'errore_relativo', 'dice3d'}), fullfile(outDir, 'volumes.csv'));

%% Figure: distribuzione, Dice lungo il tumore, volumi
f = figure('Visible', 'off', 'Position', [0 0 1100 380]); theme(f, 'light');
subplot(1, 2, 1); histogram(D3, 0:0.05:1); xlabel('volume Dice'); ylabel('patients'); title(sprintf('Test set (%d patients)', n));
subplot(1, 2, 2); plot(1:5, Spos.dice_slice_mediana, 'o-'); xticks(1:5); xticklabels(lab); ylim([0 1]); grid on;
xlabel('position along the tumour'); ylabel('slice Dice (median)'); title('Quality along the tumour');
exportgraphics(f, fullfile(outDir, 'test_distribution.png'), 'Resolution', 110); savePanels(f, fullfile(outDir, 'panels', 'test_distribution')); close(f);
f = figure('Visible', 'off', 'Position', [0 0 1000 440]); theme(f, 'light');
subplot(1, 2, 1);                                    % volume stimato contro esperto
loglog(VE, max(VP, 0.1), '.', 'MarkerSize', 8); hold on; lim = [0.1 max([VE; VP]) * 1.2]; loglog(lim, lim, 'k--');
xlabel('expert volume (cm^3)'); ylabel('estimated volume (cm^3)'); grid on; axis([lim lim]); axis square;
title(sprintf('Volumes: Spearman %.2f, median error %+.0f%%', rs, 100 * median(rel)), 'FontSize', 9);
subplot(1, 2, 2);                                    % precisione contro recall per voxel (risultati vuoti esclusi)
plot(SENS, PREC, '.', 'MarkerSize', 8); hold on; plot([0 1], [0 1], 'k--');
plot(median(SENS), median(PREC, 'omitnan'), 'rp', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
xlabel('voxel recall'); ylabel('voxel precision'); grid on; axis([0 1 0 1]); axis square;
title(sprintf('Median precision %.2f, recall %.2f (star)', median(PREC, 'omitnan'), median(SENS)), 'FontSize', 9);
exportgraphics(f, fullfile(outDir, 'volumes_precision_recall.png'), 'Resolution', 120); savePanels(f, fullfile(outDir, 'panels', 'volumes_precision_recall')); close(f);

%% Figura: caso migliore, mediano e peggiore, cinque slice ciascuno
[~, o] = sort(D3); pick = [o(end), o(round(n / 2)), o(1)]; tags = ["best", "median", "worst"];
f = figure('Visible', 'off', 'Position', [0 0 1300 800]); theme(f, 'light');
tl = tiledlayout(f, 3, 5, 'TileSpacing', 'compact');
for c = 1:3
    i = pick(c);
    [I, S, sp] = readCase(T.img(i), T.seg(i)); ks = find(squeeze(any(S, [1 2])));
    [M3, k0] = segmentCase(int16(I), S, sp);
    [ux, uy] = find(any(S, 3) | any(M3, 3)); h = round(20 / sp(1));
    x = max(min(ux) - h, 1):min(max(ux) + h, size(I, 1)); y = max(min(uy) - h, 1):min(max(uy) + h, size(I, 2));
    kk = unique(round(linspace(ks(1), ks(end), 5))); kk(abs(kk - k0) == min(abs(kk - k0))) = k0;   % include la partenza
    for k = kk
        nexttile(tl); imshow(mat2gray(I(x, y, k)'), []); hold on;
        if any(S(:, :, k), 'all'), visboundaries(S(x, y, k)', 'Color', 'r', 'LineWidth', 0.6); end
        if any(M3(:, :, k), 'all'), visboundaries(M3(x, y, k)', 'Color', 'g', 'LineWidth', 0.6); end
        id = char(ids(i)); start = ""; if k == k0, start = " (start)"; end
        title(sprintf('%s %s (Dice %.2f), slice %d%s', tags(c), id(end-5:end), D3(i), k, start), 'FontSize', 7);
    end
    for e = numel(kk) + 1:5, nexttile(tl); axis off; end
end
exportgraphics(f, fullfile(outDir, 'cases_best_median_worst.png'), 'Resolution', 110); savePanels(f, fullfile(outDir, 'panels', 'cases_best_median_worst')); close(f);
fprintf('\nFatto in %.0f s. Risultati in %s\n', toc, outDir);
end

function r = ranks(x)
% Ranghi (1 = valore piu' piccolo), al posto di tiedrank (Statistics Toolbox).
% ponytail: pareggi non mediati, rari con volumi continui; mediarli se servono ranghi esatti.
    [~, o] = sort(x); r = zeros(size(x)); r(o) = 1:numel(x);
end
