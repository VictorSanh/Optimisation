
function[F, G, H, ind] = OraclePH(qc, ind)

    // Si ind = 2, on calcule seulement F

    if ind == 2 then 

        F = 1/3*(q0+B*qc)'*(r .* (q0+B*qc) .* abs(q0+B*qc)) + pr'*Ar*(q0+B*qc) ;

        // Si ind = 3, on calcule seulement G

    elseif ind == 3 then

        G = B'*Ar'*pr + B'*(r.*(q0+B*qc).*abs(q0+B*qc));
        F = 0

        // Si ind = 4, on calcule F et G 
    elseif ind == 4 then

        F =  1/3*(q0+B*qc)'*(r .* (q0+B*qc) .* abs(q0+B*qc)) + (pr'*Ar*(q0+B*qc));
        G = B'*Ar'*pr + B'*(r .* (q0+B*qc) .* abs(q0+B*qc));

        // Si ind = 5, on calcule seulement H     
    elseif ind == 5 then

        temp = (r.*abs(q0+B*qc))*ones(1, 9);
        H = 2*B'*(temp.*B);
        F = 0
        G = 0

        // Si ind = 6, on calcule G et H 
    elseif ind == 6 then

        temp = (r.*abs(q0+B*qc))*ones(1, 9);
        H = 2*B'*(temp.*B);
        G = B'*Ar'*pr + B'*(r .* (q0+B*qc) .* abs(q0+B*qc));
        F = 0

        // Si ind = 7, on calcule F, G et H 
    elseif ind == 7 then

        temp = (r.*abs(q0+B*qc))*ones(1, size(B, 2));
        H = 2*B'*(temp.*B);
        F =  1/3*(q0+B*qc)'*(r .* (q0+B*qc) .* abs(q0+B*qc)) + (pr'*Ar*(q0+B*qc));
        G = B'*Ar'*pr + B'*(r .* (q0+B*qc) .* abs(q0+B*qc));

    end
endfunction

