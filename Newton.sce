function [fopt,xopt,gopt]=Newton(Oracle,xini)


    ///////////////////////////////////////////////////////////////////////////////
    //                                                                           //
    //         RESOLUTION D'UN PROBLEME D'OPTIMISATION SANS CONTRAINTES          //
    //                                                                           //
    //         Methode de Newton                                                 //
    //                                                                           //
    ///////////////////////////////////////////////////////////////////////////////

    exec('Wolfe_Skel.sci')
    // ------------------------
    // Parametres de la methode
    // ------------------------

    titre = "Parametres de la méthode de Newton";
    labels = ["Nombre maximal d''iterations";...
    "Valeur du pas de gradient";...
    "Seuil de convergence sur ||G||"];
    typ = list("vec",1,"vec",1,"vec",1);
    default = ["5000";"1.0";"0.0000000001"];
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
    for k = 1:iter

        //    - valeur du critere et du gradient

        ind = 7;
        [F,G,H] = Oracle(x,ind);

        //    - test de convergence

        if norm(G) <= tol then
            kstar = k;
            break
        end

        //    - calcul de la direction de descente
        //On vérifie tout de même que la matrice hessienne est inversible...
        //Si ce n'est pas le cas, on sort de la boucle...
        if ~(det(H) == 0) then
            D = -inv(H)*G
        else
            disp('Attention ! Hessienne non inversible !')
            break
        end

        //Calcul de la longueur du pas de gradient
        alphan = Wolfe(alpha, x, D, Oracle);

        //    - mise a jour des variables
        x = x + (alphan*D);
        alpha = alphan

        //    - evolution du gradient, du pas et du critere
        logG = [ logG ; log10(norm(G)) ];
        logP = [ logP ; log10(alpha) ];
        Cout = [ Cout ; F ];

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
    disp('Fin de la methode de Newton')
    disp(cvge)

    // - visualisation de la convergence

    Visualg(logG,logP,Cout);

endfunction
