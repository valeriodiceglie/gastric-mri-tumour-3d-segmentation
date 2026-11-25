function renderTestCases3D()
%RENDERTESTCASES3D Ricostruzione 3D (marching cubes) del caso migliore, mediano e peggiore del test.
%   I casi sono scelti dal Dice per volume di results/analyzeTest/per_patient.csv (eseguire prima
%   analyzeTest). Per ciascuno: risultato | esperto | sovrapposizione con il piano della slice di
%   partenza. Assi in mm (voxel non cubici). Le superfici sono lisciate (smooth3, gaussiano 3x3x3)
%   solo per la visualizzazione: Dice e volumi usano le maschere originali.
%   Uscita: results/analyzeTest/render3d_best_median_worst.png
T = readCases("manifest/test_cases.csv"); T = T(T.modality == "V", :);
outDir = "results/analyzeTest";
R = readCases(fullfile(outDir, 'per_patient.csv'));
[~, o] = sort(str2double(R.dice3d)); pick = R.patient_id([o(end), o(round(height(R) / 2)), o(1)]);
tags = ["Best", "Median", "Worst"];
gray = [0.78 0.78 0.80]; red = [0.85 0.15 0.15]; green = [0.10 0.70 0.25]; blue = [0.2 0.4 0.9];

f = figure('Visible', 'off', 'Position', [0 0 1350 1250]); theme(f, 'light');
tl = tiledlayout(f, 3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for c = 1:3
    r = T(T.patient_id == pick(c), :);
    [I, S, sp] = readCase(r.img, r.seg);
    [M3, k0] = segmentCase(int16(I), S, sp);
    % ritaglio attorno a esperto e risultato
    [ux, uy, uz] = ind2sub(size(S), find(S | M3)); m = 3;
    x = max(min(ux) - m, 1):min(max(ux) + m, size(S, 1)); y = max(min(uy) - m, 1):min(max(uy) + m, size(S, 2));
    z = max(min(uz) - 1, 1):min(max(uz) + 1, size(S, 3));
    Sc = S(x, y, z); Mc = M3(x, y, z);
    id = char(pick(c));
    hdr = sprintf('%s %s: Dice %.2f, expert %.1f cm^3, estimated %.1f cm^3', tags(c), id(end-5:end), ...
        diceCoef(M3, S), nnz(S) * prod(sp) / 1000, nnz(M3) * prod(sp) / 1000);
    ax = nexttile(tl); drawSurface(ax, Mc, sp, green, 1); title(ax, {hdr, 'result'}, 'FontSize', 8);
    ax = nexttile(tl); drawSurface(ax, Sc, sp, gray, 1); title(ax, 'expert', 'FontSize', 8);
    ax = nexttile(tl); drawSurface(ax, Sc, sp, red, 0.25); drawSurface(ax, Mc, sp, green, 0.85);
    [ny, nx] = size(Sc, [1 2]); zp = (k0 - z(1) + 1) * sp(3);
    patch(ax, [0 nx nx 0] * sp(2), [0 0 ny ny] * sp(1), zp * [1 1 1 1], blue, 'FaceAlpha', 0.12, 'EdgeColor', blue);
    title(ax, 'overlay (red expert, green result, blue = start slice)', 'FontSize', 8);
end
exportgraphics(f, fullfile(outDir, 'render3d_best_median_worst.png'), 'Resolution', 130); savePanels(f, fullfile(outDir, 'panels', 'render3d_best_median_worst')); close(f);
fprintf('Figura in %s\n', outDir);
end

function drawSurface(ax, M, sp, col, alpha)
% Superficie con marching cubes (isosurface), aggiunta agli assi. Vertici in mm: x = colonne, y = righe, z = slice.
    hold(ax, 'on');
    if ~any(M(:))
        text(ax, 0.5, 0.5, 0.5, 'no segmentation', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 9);
    else
        fv = isosurface(smooth3(double(padarray(M, [1 1 1])), 'gaussian', 3), 0.5);
        fv.vertices = (fv.vertices - 1) .* [sp(2) sp(1) sp(3)];
        patch(ax, fv, 'FaceColor', col, 'EdgeColor', 'none', 'FaceAlpha', alpha);
    end
    daspect(ax, [1 1 1]); view(ax, -37.5, 25); axis(ax, 'tight'); grid(ax, 'on'); box(ax, 'on');
    xlabel(ax, 'x (mm)'); ylabel(ax, 'y (mm)'); zlabel(ax, 'z (mm)');
    camlight(ax, 'headlight'); lighting(ax, 'gouraud'); material(ax, 'dull');
end
