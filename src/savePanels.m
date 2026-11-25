function savePanels(f, stem)
%SAVEPANELS Salva ogni riquadro (asse) della figura f in un PNG a se', per le sottofigure della relazione.
%   File: <stem>_1.png, <stem>_2.png, ... nell'ordine in cui i riquadri sono stati creati.
%   I titoli dei riquadri vengono tolti: nella relazione il loro contenuto va nelle didascalie.
%   Va chiamata dopo aver salvato la figura completa (che tiene i titoli). La cartella viene creata se manca.
    [~, ~] = mkdir(fileparts(stem));
    ax = flip(findobj(f, 'Type', 'axes'));
    for i = 1:numel(ax)
        title(ax(i), '');
        exportgraphics(ax(i), sprintf('%s_%d.png', stem, i), 'Resolution', 200);
    end
end
