function[F, G, H, ind] = OracleDH(lambda, ind)

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

        // Si ind = 5, on calcule seulement H     
    elseif ind == 5 then
        
 H = 1/2*Ad./sqrt(r*ones(1,size(Ad,1)).*(abs(zH)*ones(1,size(Ad,1))))' * Ad' ;
        //H = 2*(r.*abs(q)*ones(1,size(r,1))*eye(size(r,1), size(r,1)));
        F = 0
        G = 0

        // Si ind = 6, on calcule G et H 
    elseif ind == 6 then

H = 1/2*Ad./sqrt(r*ones(1,size(Ad,1)).*(abs(zH)*ones(1,size(Ad,1))))' * Ad' ;
        G = -(Ad*qH-fd);
        F = 0

        // Si ind = 7, on calcule F, G et H 
    elseif ind == 7 then

H = 1/2*Ad./sqrt(r*ones(1,size(Ad,1)).*(abs(zH)*ones(1,size(Ad,1))))' * Ad' ;
        //H = -2*(r.*abs(q)*ones(1,size(r,1))*eye(size(r,1), size(r,1)));
        F = -(1/3 * (zH' * qH) + (pr' * Ar * qH) + (lambda' * (Ad*qH - fd)));
        G = -(Ad*qH-fd);

    end
endfunction

