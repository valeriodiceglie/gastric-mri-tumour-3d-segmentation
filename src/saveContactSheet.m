function saveContactSheet(file, tiles, header)
%SAVECONTACTSHEET Figura con un riquadro per paziente: slice, contorno esperto e risultato.
%   tiles(i) : img (slice 2D), ref (maschera esperta), mask (risultato), box [x0 x1 y0 y1] (ritaglio),
%              seeds (N x 2), label (titolo del riquadro)
%   Rosso = esperto, verde = risultato, + giallo = semi. Slice in vista radiologica (trasposta).
    f = figure('Visible', 'off', 'Position', [0 0 1300 1100]); theme(f, 'light');
    tl = tiledlayout(f, 4, ceil(numel(tiles) / 4), 'TileSpacing', 'compact');
    for t = tiles(:)'
        b = t.box; r = b(1):b(2); c = b(3):b(4);
        nexttile(tl); imshow(mat2gray(t.img(r, c)'), []); hold on;
        visboundaries(t.ref(r, c)', 'Color', 'r', 'LineWidth', 0.6);
        if any(t.mask(r, c), 'all'), visboundaries(t.mask(r, c)', 'Color', 'g', 'LineWidth', 0.6); end
        plot(t.seeds(:, 1) - b(1) + 1, t.seeds(:, 2) - b(3) + 1, 'y+', 'MarkerSize', 6);
        title(t.label, 'FontSize', 8);
    end
    title(tl, header + " (red expert, green result, + seeds)", 'FontSize', 9, 'Interpreter', 'none');
    exportgraphics(f, file, 'Resolution', 100); close(f);
end
