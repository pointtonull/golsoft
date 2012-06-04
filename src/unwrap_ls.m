function unwrapped=unwrap_ls(wrapped)
    [rows, cols] = size(wrapped);

    for row=1:rows
        for col=1:cols
            if row < rows
                difefila(row, col) = wrapped(row + 1, col) - wrapped(row, col);
            else
                difefila(row, col) = 0;
            end
        end
    end

    for row = 1:rows
        for col = 1:cols
            if col < cols
                difecol(row, col) = wrapped(row, col + 1) - wrapped(row, col);
            else
                difecol(row, col) = 0;
            end
        end
    end

    for row = 1:rows
        for col = 1:cols
            W(row, col) = atan2(sin(difefila(row, col)),
                cos(difefila(row, col)));
        end
    end

    for row = 1:rows
        for col = 1:cols
            V(row, col) = atan2(sin(difecol(row, col)), cos(difecol(row, col)));
        end
    end

    for row = 1:rows
        for col = 1:cols
            if row < 2
                rho_x(row, col) = W(row, col);
            else
                rho_x(row, col)= W(row, col) - W(row - 1, col);
            end
         end
    end

    for row = 1:rows
        for col = 1:cols
            if col < 2
                rho_y(row, col) = V(row, col);
            else
                rho_y(row, col) = V(row, col) - V(row, col - 1);
            end
         end
    end

    rho = rho_x + rho_y;
    DCT = dct2(rho);

    for row = 1:rows
        for col = 1:cols
            cosines(row, col) = (2 * (cos(pi * row / rows)
                + cos(pi * col / cols) - 2));
        end
    end
    phiINV = DCT ./ cosines;
    unwrapped = idct2(phiINV);
end
