///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  MONITEUR D'ENCHAINEMENT POUR LE CALCUL DE L'EQUILIBRE D'UN RESEAU D'EAU  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

// --------------------------------------
// Dimensionnement de l'espace de travail
// --------------------------------------

stacksize(10000000);

// ------------------------------------------
// Fonctions fournies dans le cadre du projet
// ------------------------------------------

// Donnees du problemes

exec('Probleme_R.sce');
exec('Structures_R.sce');

// Affichage des resultats

exec('Visualg.sci');

// Verification  des resultats

exec('HydrauliqueP.sci');
exec('HydrauliqueD.sci');
exec('Verification.sci');

// ------------------------------------------
// Fonctions a ecrire dans le cadre du projet
// ------------------------------------------

// ---> Charger les fonctions  associees a l'oracle du probleme,
//      aux algorithmes d'optimisation et de recherche lineaire.
//
// Exemple : la fonction "optim" de Scilab
//
exec('OraclePG.sci');
exec('OraclePH.sce');
exec('OracleDG.sce');
exec('OracleDH.sce');

exec('Gradient_F.sci')
exec('Gradient_V.sce');
exec('Gradient_Conj.sce');
exec('Grad_QN.sce');
exec('Newton.sce');

exec('Optim_Scilab.sci');

titrgr = "Gradient à pas variable sur le problème dual"
//titrgr = "Fonction optim de Scilab sur le probleme primal";
//titrgr = "Fonction optim de Scilab sur le probleme dual";

// ------------------------------
// Initialisation de l'algorithme
// ------------------------------

// La dimension (n-md) est celle du probleme primal

xini = 0.1 * rand(n-md,1);
lambda = 0.1*rand(md,1);

// ----------------------------
// Minimisation proprement dite
// ----------------------------

// Exemple : la fonction "optim" de Scilab
//[fopt,xopt,gopt] = Optim_Scilab(OraclePG,xini);

//[fopt, xopt, gopt] = Gradient_F(OraclePG, xini);
//[fopt, xopt, gopt] = Gradient_F(OracleDG, lambda);
//[fopt, xopt, gopt] = Gradient_V(OraclePG, xini);
//[fopt, xopt, gopt] = Gradient_V(OracleDG, lambda);
//[fopt, xopt, gopt] = Gradient_Conj(OraclePG, xini);
//[fopt, xopt, gopt] = Gradient_Conj(OracleDG, lambda);
//[fopt, xopt, gopt] = Gradient_QN(OraclePG, xini);
//[fopt, xopt, gopt] = Gradient_QN(OracleDG, lambda);

//[fopt, xopt, gopt] = Newton(OraclePH, xini);
//[fopt, xopt, gopt] = Newton(OracleDH, lambda);

// --------------------------
// Verification des resultats
// --------------------------

[q,z,f,p] = HydrauliqueP(xopt);
//[q,z,f,p] = HydrauliqueD(xopt);

Verification(q,z,f,p);

//
