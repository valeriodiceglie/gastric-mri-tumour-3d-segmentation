function plotPipelineSteps()
%PLOTPIPELINESTEPS Figura dei passi del metodo finale su un paziente dev tipico (per la relazione).
%   Paziente = quello con il Dice 2D mediano sulla slice di partenza (nessuna scelta a mano).
%   Riga 1: i passi di segmentSlice sulla slice di partenza (ROI e semi, normalizzazione e filtro,
%   soglia adattiva, region growing, apertura e riempimento). Riga 2: cinque slice consecutive
%   attorno alla partenza dopo la propagazione (rosso esperto, verde risultato).
%   Seconda figura: un passo di propagazione, dalla slice di partenza k0 a k0+1 (maschera
%   precedente, marker eroso, segmentazione dal marker, vincolo di 1 mm, maschera accettata).
%   Uscite: results/plotPipelineSteps/pipeline_steps.png, propagation_step.png
T = readCases('manifest/dev_cases.csv');
outDir = 'results/plotPipelineSteps'; [~, ~] = mkdir(outDir);
[P, Q] = methodConfig();

%% Paziente con il Dice 2D mediano sulla slice di partenza
C = loadStartSlices(T, "V"); d = zeros(numel(C), 1);
for i = 1:numel(C)
    d(i) = diceCoef(segmentSlice(C(i).img, pickSeeds(C(i).ref, C(i).sp, Q.minSeedAreaMM2), C(i).sp, P), C(i).ref);
end
[~, o] = sort(d); c = C(o(ceil(numel(C) / 2)));
seeds = pickSeeds(c.ref, c.sp, Q.minSeedAreaMM2);
[mask, box, steps] = segmentSlice(c.img, seeds, c.sp, P);
r = box(1):box(2); q = box(3):box(4); s = seeds - [box(1) - 1, box(3) - 1];

f = figure('Visible', 'off', 'Position', [0 0 1500 640]); theme(f, 'light');
tl = tiledlayout(f, 2, 5, 'TileSpacing', 'compact', 'Padding', 'compact');

%% Riga 1: passi sulla slice di partenza (ROI, in vista radiologica)
m = round(numel(r) / 2);                                      % margine attorno alla ROI nel primo pannello
rr = max(box(1) - m, 1):min(box(2) + m, size(c.img, 1)); qq = max(box(3) - m, 1):min(box(4) + m, size(c.img, 2));
nexttile(tl); imshow(mat2gray(c.img(rr, qq)'), []); hold on;
rectangle('Position', [box(1) - rr(1) + 0.5, box(3) - qq(1) + 0.5, numel(r), numel(q)], 'EdgeColor', 'y', 'LineWidth', 1);
plot(seeds(:, 1) - rr(1) + 1, seeds(:, 2) - qq(1) + 1, 'y+', 'MarkerSize', 9, 'LineWidth', 1.5);
title('start slice, 60 mm ROI, seeds');
panel(tl, steps.roi, s, 'normalised + Gaussian');
panel(tl, steps.threshold, s, 'adaptive threshold');
panel(tl, steps.grown, s, 'region growing');
nexttile(tl); imshow(mat2gray(c.img(r, q)'), []); hold on;
visboundaries(c.ref(r, q)', 'Color', 'r', 'LineWidth', 0.8); visboundaries(mask(r, q)', 'Color', 'g', 'LineWidth', 0.8);
title(sprintf('opening + filling, Dice %.2f', diceCoef(mask, c.ref)));

%% Riga 2: propagazione, cinque slice consecutive attorno alla partenza
r3 = T(T.patient_id == c.id & T.modality == "V", :);
[I, S, sp] = readCase(r3.img, r3.seg);
M3 = segmentCase(int16(I), S, sp);
for k = c.k + (-2:2)
    nexttile(tl); imshow(mat2gray(I(r, q, k)'), []); hold on;
    if any(S(r, q, k), 'all'), visboundaries(S(r, q, k)', 'Color', 'r', 'LineWidth', 0.8); end
    if any(M3(r, q, k), 'all'), visboundaries(M3(r, q, k)', 'Color', 'g', 'LineWidth', 0.8); end
    title(sprintf('slice k_0%+d, Dice %.2f', k - c.k, diceCoef(M3(:, :, k), S(:, :, k))));
end
exportgraphics(f, fullfile(outDir, 'pipeline_steps.png'), 'Resolution', 200); savePanels(f, fullfile(outDir, 'panels', 'pipeline_steps')); close(f);

%% Seconda figura: un passo di propagazione k0 -> k0+1, con i parametri di propagateVolume
% ponytail: ricostruisce a mano i primi passi di propagateVolume per mostrarli; il pannello (e)
% e' invece il risultato vero di segmentCase, quindi un'eventuale differenza resta visibile.
k1 = c.k + 1; pxmm = mean(sp(1:2)); prev = M3(:, :, c.k); img = double(I(:, :, k1));
marker = imerode(prev, strel('disk', round(Q.erodeMM / pxmm), 0));
Pk = P; Pk.perSeed = false;
cand = segmentSlice(img, marker, sp(1:2), Pk);
band = imdilate(prev, strel('disk', round(Q.clipMM / pxmm), 0));
f = figure('Visible', 'off', 'Position', [0 0 1500 340]); theme(f, 'light');
tl = tiledlayout(f, 1, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
show = @() imshow(mat2gray(img(r, q)'), []);
nexttile(tl); show(); hold on; visboundaries(prev(r, q)', 'Color', 'y', 'LineWidth', 0.8);
title('slice k_0+1, previous mask M_{k_0}');
nexttile(tl); show(); hold on; visboundaries(marker(r, q)', 'Color', 'c', 'LineWidth', 0.8);
title('marker: M_{k_0} eroded by 1 mm');
nexttile(tl); show(); hold on; visboundaries(cand(r, q)', 'Color', 'm', 'LineWidth', 0.8);
title('same steps from the marker');
nexttile(tl); show(); hold on; visboundaries(band(r, q)', 'Color', 'y', 'LineStyle', '--', 'LineWidth', 0.8);
visboundaries(cand(r, q)' & band(r, q)', 'Color', 'g', 'LineWidth', 0.8);
title('kept within 1 mm of M_{k_0}');
nexttile(tl); show(); hold on;
visboundaries(S(r, q, k1)', 'Color', 'r', 'LineWidth', 0.8); visboundaries(M3(r, q, k1)', 'Color', 'g', 'LineWidth', 0.8);
title(sprintf('accepted M_{k_0+1}, Dice %.2f', diceCoef(M3(:, :, k1), S(:, :, k1))));
exportgraphics(f, fullfile(outDir, 'propagation_step.png'), 'Resolution', 200); savePanels(f, fullfile(outDir, 'panels', 'propagation_step')); close(f);
fprintf('Passo k0->k0+1: area candidato %d, dopo il vincolo %d, accettata %d pixel\n', nnz(cand), nnz(cand & band), nnz(M3(:, :, k1)));
fprintf('Paziente %s, slice %d, Dice 2D %.3f, Dice per volume %.3f\n', c.id, c.k, diceCoef(mask, c.ref), diceCoef(M3, S));
end

function panel(tl, img, seeds, label)
% Un pannello della ROI con i semi.
    nexttile(tl); imshow(img', []); hold on;
    plot(seeds(:, 1), seeds(:, 2), 'y+', 'MarkerSize', 9, 'LineWidth', 1.5);
    title(label);
end
