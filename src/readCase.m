function [I, S, sp] = readCase(imgPath, segPath)
%READCASE Legge volume e maschera NIfTI (.nii.gz) di un caso GCM 2025.
%   I  : volume double [X Y Z];  S : maschera logica [X Y Z];  sp : spaziatura [sx sy sz] in mm.
%   Orientamento LPS: per visualizzare una slice in vista radiologica usare I(:,:,k)'.
    imgPath = char(imgPath); segPath = char(segPath);
    I    = double(niftiread(imgPath));
    S    = niftiread(segPath) > 0;
    info = niftiinfo(imgPath);
    sp   = double(info.PixelDimensions(1:3));
    assert(isequal(size(I), size(S)), 'Volume e maschera hanno griglie diverse: %s', imgPath);
end
