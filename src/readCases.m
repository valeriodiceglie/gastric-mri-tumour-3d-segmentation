function T = readCases(csvPath)
%READCASES Legge manifest/<split>_cases.csv con separatore virgola e tutte le colonne come testo
%   (gli ID paziente hanno zeri iniziali che una lettura numerica perderebbe).
    opts = detectImportOptions(csvPath, 'Delimiter', ',', 'FileType', 'text');
    opts = setvartype(opts, 'string');
    T = readtable(csvPath, opts);
end
