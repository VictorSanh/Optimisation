function [fopt,xopt,gopt]=Gradient_Conj(Oracle,xini)


    ///////////////////////////////////////////////////////////////////////////////
    //                                                                           //
    //         RESOLUTION D'UN PROBLEME D'OPTIMISATION SANS CONTRAINTES          //
    //                                                                           //
    //         Methode de gradient conjugué                                      //
    //                                                                           //
    ///////////////////////////////////////////////////////////////////////////////

    exec('Wolfe_Skel.sci')
    // ------------------------
    // Parametres de la methode
    // ------------------------

    titre = "Parametres du gradient Conjugué";
    labels = ["Nombre maximal d''iterations";...
    "Valeur du pas de gradient";...
    "Seuil de convergence sur ||G||"];
    typ = list("vec",1,"vec",1,"vec",1);
    default = ["1000";"1";"0.000001"];
    [ok,iter,alpha,tol] = getvalue(titre,labels,typ,default);

    // ----------------------------
    // Initialisation des variables
    // ----------------------------

    logG = [];
    logP = [];
    Cout = [];

    timer();

    // -------------------------
    // Boucle sur les iterations
    // -------------------------

    x = xini;
    kstar = iter;
    [F,G] = Oracle(x,4);
    for k = 1:iter


        //    - valeur du critere et du gradient

        ind = 4;


        if(k == 1) then
            D = -G;
        else

            //    - test de convergence

            if norm(G) <= tol then
                kstar = k;
                break
            end

            //    - calcul de la direction de descente

            [Fn,Gn] = Oracle(x, 4);    
            beta = (Gn-G)'*Gn/(G'*G)    
            D = -Gn + beta*D

            //    - calcul de la longueur du pas de gradient

            alphan = Wolfe(alpha, x, D, Oracle);


            //    - mise a jour des variables

            x = x + (alphan*D);
            //alpha = alphan
            F = Fn
            G = Gn
            //    - evolution du gradient, du pas et du critere

            logG = [ logG ; log10(norm(G)) ];
            logP = [ logP ; log10(alphan) ];
            Cout = [ Cout ; F ];

        end

    end

    // ---------------------------
    // Resultats de l'optimisation
    // ---------------------------

    fopt = F;
    xopt = x;
    gopt = G;

    tcpu = timer();

    cvge = ['Iteration         : ' string(kstar);...
    'Temps CPU         : ' string(tcpu);...
    'Critere optimal   : ' string(fopt);...
    'Norme du gradient : ' string(norm(gopt))];
    disp('Fin de la methode de gradient conjugué')
    disp(cvge)

    // - visualisation de la convergence

    Visualg(logG,logP,Cout);

endfunction
