function unwrapped=unwrap_wls(wrapped)
    tau = 2 * pi;
    [rows, cols] = size(wrapped);

    rowdiff = [diff(wrapped, 1, 1); zeros(1, cols)];
    coldiff = [diff(wrapped, 1, 2), zeros(rows, 1)];
    
    wrowdiff = [mod(rowdiff + pi, tau)] - pi;
    wcoldiff = [mod(coldiff + pi, tau)] - pi;
    
    rhox = diff([zeros(1, rows); wrowdiff], 1, 1);
    rhoy = diff([zeros(cols, 1), wcoldiff], 1, 2);
    
    rho = rhox + rhoy;
    dct = dct2(rho);

    row = repmat((1: rows)' / rows, 1, rows);
    col = repmat((1: cols)  / cols, cols, 1);
    cosines = (2 * (cos(pi * row) + cos(pi * col) - 2));

    phiinv = dct ./ cosines;
    unwrapped = idct2(phiinv);
end
