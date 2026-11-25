function txt = describeParams(s, names)
%DESCRIBEPARAMS Testo "campo=valore ..." dei campi scelti di una struct, per tabelle e titoli.
    parts = strings(1, numel(names));
    for f = 1:numel(names), parts(f) = names{f} + "=" + string(s.(names{f})); end
    txt = strjoin(parts, ' ');
end
