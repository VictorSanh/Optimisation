
function[F, G, H, ind] = OraclePH(qc, lambda, ind)

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

        // Si ind = 5, on calcule seulement H     
    elseif ind == 5 then

        H = 2*(r.*abs(q)*ones(1,size(r,1))*eye(size(r,1), size(r,1));
        F = 0
        G = 0

        // Si ind = 6, on calcule G et H 
    elseif ind == 6 then

        H = 2*(r.*abs(q)*ones(1,size(r,1))*eye(size(r,1), size(r,1));
        G = r .* q .* abs(q) - Ar'*pr - lambda*Ad;
        F = 0

        // Si ind = 7, on calcule F, G et H 
    elseif ind == 7 then

        H = 2*(r.*abs(q)*ones(1,size(r,1))*eye(size(r,1), size(r,1));
        F = 1/3*q'*(r .* q .* abs(q)) - pr'*Ar*q - lambda*Ad*q + lambda*fd;
        G = r .* q .* abs(q) - Ar'*pr - lambda*Ad;

    end
endfunction

