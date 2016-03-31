function[F, G, ind] = OraclePG(qc, ind)

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

    end
endfunction

