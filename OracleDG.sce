function[F, G, ind] = OracleDG(lambda, ind)

    zH = -(Ar'*pr+Ad'*lambda)
    qH = sign(zH).*sqrt( abs(zH./r) ) 

    // Si ind = 2, on calcule seulement F
    if ind == 2 then 
        F = -(1/3 * (zH' * qH) + (pr' * Ar * qH) + (lambda' * (Ad*qH - fd)));

    // Si ind = 3, on calcule seulement G

    elseif ind == 3 then

        G = -(Ad*qH-fd);
        F = 0

    // Si ind = 4, on calcule F et G 

    elseif ind == 4 then
        F = -(1/3 * (zH' * qH) + (pr' * Ar * qH) + (lambda' * (Ad*qH - fd)));
        G = -(Ad*qH-fd);

    end
endfunction

