function ablationRobustness(split)
%ABLATIONROBUSTNESS Contributo di ogni passo del metodo, robustezza e casi di errore (modalita' V).
%   split : "dev" o "test". Nessun parametro viene tarato qui: tutto viene da methodConfig.
%   1) Ablazioni: i passi si aggiungono uno alla volta e si misura il Dice per volume:
%        1 Otsu, un seme, solo slice di partenza   2 + soglia adattiva   3 + morfologia
%        4 + piu' semi   5 + propagazione senza vincoli   6 + propagazione con vincoli (metodo finale)
%      con IC 95% della mediana (bootstrap) e Wilcoxon appaiato contro il passo precedente.
%   2) Robustezza: slice di partenza spostata di -5..+5 (semi dall'esperto su quella slice);
%      semi spostati di 2, 4, 6 mm nelle quattro direzioni (media sulle direzioni).
%   3) Casi di errore: slice segmentate dove l'esperto non vede tumore.
%   Un paziente alla volta, per contenere la memoria.
%   Uscite in results/ablationRobustness/<split>/: ablation_summary.csv, ablation_per_patient.csv,
%   robustness_start_slice.csv, robustness_seed_shift.csv, ablation_robustness.png, error_cases.png
if nargin < 1, split = "dev"; end
T = readCases("manifest/" + split + "_cases.csv");
T = sortrows(T(T.modality == "V", :), 'patient_id'); n = height(T);
outDir = fullfile("results", "ablationRobustness", split); [~, ~] = mkdir(outDir); tic;

%% Varianti delle ablazioni
[P_final, Q_final, Q_v1] = methodConfig();
P_adapt = P_final; P_adapt.openMM = 0; P_adapt.fill = false;            % senza morfologia
P_otsu = P_adapt; P_otsu.method = 'otsu'; P_otsu.sigmaMM = 1; P_otsu.conn = 8;
variants = ["1 Otsu, un seme, solo slice di partenza"; "2 + soglia adattiva"; "3 + morfologia"; ...
            "4 + piu' semi"; "5 + propagazione senza vincoli"; "6 + propagazione con vincoli (metodo finale)"];
OFFSETS = -5:5; SHIFTS_MM = [2 4 6]; DIRS = [1 0; -1 0; 0 1; 0 -1];
A = nan(n, numel(variants)); A2d = nan(n, 4); ROB = nan(n, numel(OFFSETS)); SH = nan(n, numel(SHIFTS_MM));
nFP = zeros(n, 1); nTumor = zeros(n, 1);
err = struct('id', {}, 'img', {}, 'pred', {}, 'dist', {}, 'nfp', {});

for i = 1:n
    [I, S, sp] = readCase(T.img(i), T.seg(i)); I = int16(I); sp2 = sp(1:2);
    area = squeeze(sum(S, [1 2])); [~, k0] = max(area); ks = find(area); nTumor(i) = numel(ks);
    ref0 = S(:, :, k0); seed1 = centroidSeed(ref0); seedsAll = pickSeeds(ref0, sp2, Q_final.minSeedAreaMM2);
    slice0 = double(I(:, :, k0));

    %% 1) Ablazioni
    m2 = {segmentSlice(slice0, seed1, sp2, P_otsu), segmentSlice(slice0, seed1, sp2, P_adapt), ...
          segmentSlice(slice0, seed1, sp2, P_final), segmentSlice(slice0, seedsAll, sp2, P_final)};
    for v = 1:4, A(i, v) = diceCoef(onSlice(m2{v}, size(S), k0), S); A2d(i, v) = diceCoef(m2{v}, ref0); end
    A(i, 5) = diceCoef(propagateVolume(I, sp, k0, seedsAll, P_final, Q_v1), S);
    M3 = propagateVolume(I, sp, k0, seedsAll, P_final, Q_final); A(i, 6) = diceCoef(M3, S);

    %% 3) Casi di errore: la slice senza tumore con l'area segmentata piu' grande
    fp = find(squeeze(any(M3, [1 2])) & area == 0); nFP(i) = numel(fp);
    if ~isempty(fp)
        [~, w] = max(squeeze(sum(M3(:, :, fp), [1 2]))); kf = fp(w);
        [px, py] = find(M3(:, :, kf)); h = round(25 / sp(1));
        x = max(min(px) - h, 1):min(max(px) + h, size(I, 1)); y = max(min(py) - h, 1):min(max(py) + h, size(I, 2));
        err(end+1) = struct('id', T.patient_id(i), 'img', double(I(x, y, kf)), 'pred', M3(x, y, kf), ...
                            'dist', min(abs(ks - kf)), 'nfp', numel(fp)); %#ok<AGROW>
    end

    %% 2) Robustezza: slice di partenza spostata
    for o = 1:numel(OFFSETS)
        k = k0 + OFFSETS(o);
        if k < 1 || k > size(S, 3) || area(k) == 0, continue; end
        sd = pickSeeds(S(:, :, k), sp2, Q_final.minSeedAreaMM2);
        ROB(i, o) = diceCoef(propagateVolume(I, sp, k, sd, P_final, Q_final), S);
    end
    %% 2) Robustezza: semi spostati
    for s = 1:numel(SHIFTS_MM)
        d = zeros(size(DIRS, 1), 1);
        for q = 1:size(DIRS, 1)
            sd = min(max(seedsAll + round(SHIFTS_MM(s) * DIRS(q, :) ./ sp2), 1), [size(I, 1) size(I, 2)]);
            d(q) = diceCoef(propagateVolume(I, sp, k0, sd, P_final, Q_final), S);
        end
        SH(i, s) = mean(d);
    end
    fprintf('%s %d/%d  Dice3D finale %.3f  (%.0f s)\n', T.patient_id(i), i, n, A(i, 6), toc);
end

%% Riepilogo delle ablazioni: bootstrap e Wilcoxon appaiato (ogni passo contro il precedente)
ci = zeros(numel(variants), 2); pw = nan(numel(variants), 1);
for v = 1:numel(variants)
    ci(v, :) = bootMedianCI(A(:, v));
    if v > 1, pw(v) = wilcoxonSignedRank(A(:, v) - A(:, v - 1)); end
end
Sab = table(variants, median(A)', ci(:, 1), ci(:, 2), mean(A)', mean(A < 0.1)', pw, 'VariableNames', ...
    {'variante', 'dice3d_median', 'ci95_basso', 'ci95_alto', 'dice3d_mean', 'quota_fallimenti', 'p_wilcoxon_vs_precedente'});
writetable(Sab, fullfile(outDir, 'ablation_summary.csv'));
writetable([table(T.patient_id, 'VariableNames', {'patient_id'}), array2table([A A2d nTumor nFP], 'VariableNames', ...
    cellstr(["v" + (1:numel(variants)), "dice2d_start_v" + (1:4), "slice_tumore", "slice_fp_finale"]))], ...
    fullfile(outDir, 'ablation_per_patient.csv'));
fprintf('\n== Ablazioni (%s, %d pazienti) ==\n', split, n); disp(Sab);
fprintf('Dice 2D sulla slice di partenza, mediana: %s\n', mat2str(round(median(A2d), 3)));

%% Riepilogo della robustezza
Srob = table(OFFSETS', median(ROB, 'omitnan')', sum(~isnan(ROB))', 'VariableNames', ...
    {'spostamento_slice', 'dice3d_median', 'n_pazienti'});
Ssh = table(SHIFTS_MM', median(SH)', mean(SH)', 'VariableNames', {'spostamento_seme_mm', 'dice3d_median', 'dice3d_mean'});
writetable(Srob, fullfile(outDir, 'robustness_start_slice.csv'));
writetable(Ssh, fullfile(outDir, 'robustness_seed_shift.csv'));
fprintf('\n== Robustezza alla slice di partenza ==\n'); disp(Srob);
fprintf('== Robustezza allo spostamento dei semi (0 mm: %.3f) ==\n', median(A(:, 6))); disp(Ssh);
fprintf('Slice senza tumore segmentate (finale): media %.2f, pazienti con almeno una: %d su %d\n', mean(nFP), nnz(nFP), n);

%% Figure
f = figure('Visible', 'off', 'Position', [0 0 900 380]); theme(f, 'light');
subplot(1, 2, 1); boxchart(A); xticklabels(["Otsu", "+ adaptive", "+ opening", "+ seeds", "+ unconstr.", "+ constr."]);
xtickangle(30); ylabel('volume Dice'); title('Ablation');
subplot(1, 2, 2); hm = plot(OFFSETS, median(ROB, 'omitnan'), 'o-'); hold on;
hq = plot(OFFSETS, prctile(ROB, 25), ':', 'Color', [0.85 0.33 0.10]); plot(OFFSETS, prctile(ROB, 75), ':', 'Color', [0.85 0.33 0.10]); grid on;
legend([hm hq], {'median', 'quartiles (25th, 75th)'}, 'Location', 'northwest');
xlabel('start slice offset from the largest-area slice'); ylabel('volume Dice'); title('Start-slice robustness (median, quartiles)');
exportgraphics(f, fullfile(outDir, 'ablation_robustness.png'), 'Resolution', 110); savePanels(f, fullfile(outDir, 'panels', 'ablation_robustness')); close(f);
if ~isempty(err)                          % fino a 8 pazienti, quelli con piu' slice senza tumore
    [~, o] = sort([err.nfp], 'descend'); o = o(1:min(8, numel(o)));
    f = figure('Visible', 'off', 'Position', [0 0 1300 650]); theme(f, 'light');
    tl = tiledlayout(f, 2, 4, 'TileSpacing', 'compact');
    for e = o
        nexttile(tl); imshow(mat2gray(err(e).img'), []); hold on; visboundaries(err(e).pred', 'Color', 'g', 'LineWidth', 0.7);
        id = char(err(e).id); title(sprintf('%s: %d FP slices, this one %d from tumour', id(end-5:end), err(e).nfp, err(e).dist), 'FontSize', 8);
    end
    exportgraphics(f, fullfile(outDir, 'error_cases.png'), 'Resolution', 110); savePanels(f, fullfile(outDir, 'panels', 'error_cases')); close(f);
end
fprintf('\nFatto in %.0f s. Risultati in %s\n', toc, outDir);
end

function M3 = onSlice(m2, sz, k)
% Volume vuoto con la sola slice k: per il Dice per volume della sola slice di partenza.
    M3 = false(sz); M3(:, :, k) = m2;
end
