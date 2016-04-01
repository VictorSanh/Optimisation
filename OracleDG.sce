function[F, G, ind] = OracleDG(q, lambda, ind)

    // Si ind = 2, on calcule seulement F

    if ind == 2 then 

        F = 1/3*q'*(r .* q .* abs(q)) - pr'*Ar*q - lambda*Ad*q + lambda*fd;

        // Si ind = 3, on calcule seulement G

    elseif ind == 3 then

        G = r .* q .* abs(q) - Ar'*pr - lambda*Ad;
        F = 0
        // Si ind = 4, on calcule F et G 

    elseif ind == 4 then

        F = 1/3*q'*(r .* q .* abs(q)) - pr'*Ar*q - lambda*Ad*q + lambda*fd;
        G = r .* q .* abs(q) - Ar'*pr - lambda*Ad;

    end
endfunction

