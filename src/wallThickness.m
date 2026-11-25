function t = wallThickness(m, pxmm, stat)
%WALLTHICKNESS Spessore della regione m (mm) con la trasformata di distanza.
%   d = distanza di ogni pixel interno dal bordo; il diametro del cerchio inscritto in un punto e' 2*d.
%   stat = 'max'  : diametro del piu' grande cerchio inscritto (parte piu' spessa della regione)
%          'skel' : mediana di 2*d sullo scheletro (spessore tipico lungo la parete)
%   pxmm : lato del pixel (mm). Regione vuota: t = 0.
    if ~any(m(:)), t = 0; return; end
    d = bwdist(~m);
    switch stat
        case 'max'
            t = 2 * max(d(:)) * pxmm;
        case 'skel'
            sk = bwskel(m);
            if any(sk(:)), t = 2 * median(d(sk)) * pxmm; else, t = 2 * max(d(:)) * pxmm; end
    end
end
