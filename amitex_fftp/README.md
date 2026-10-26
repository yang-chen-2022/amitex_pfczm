Programme AMITEX FFTP
=====================

1 - Installation du programme
-----------------------------
Renseigner dans la variable FC le compilateur désiré (ifort ou gfortran)
Préciser dans quel dossier se trouve la librairie FFT (rien à faire pour maldives)
Vérifier au préalable que les modules nécessaires (gnu ou intel, mpi et tfel) sont chargés.
Lancer l'installation avec la commande *./install*

2 - Exécution du programme
--------------------------
Dans le dossier *cas_tests*, le script *script.sh* permet de lancer plusieurs cas tests.
Pour le lancer sur le cluster Maldives, il suffit d'exécuter la commande

    qsub script_maldives  
	
On pourra préciser dans ce script le nombre de noeuds désirés ou le noeud choisi.
Un script compatible avec le cluster Poincare de la Maison de la Simulation
est également disponible. Il faut alors lancer sur une frontale du cluster la commande
llsubmit script_poincare

3 - Graphique à partir des valeurs moyennes
-------------------------------------------
Les sorties standards (std, mstd et zstd) présentent des données moyennes
(respectivement sur la cellule, par matériau et par zone). Ces données sont 
écrites en colonnes et la grandeur correspondant à chaque colonne est précisée
en commentaire au début du fichier. Ainsi, les résultats peuvent être facilement
visualisés en utilisant gnuplot. Un exemple de script gnuplot est présenté dans 
le dossier *post/plot*. Ce dernier peut être exécuté une fois que les tests de validation
ont été lancés.

4 - Visualisation de la déformée
--------------------------------
Dans le dossier *post/deformed_shape*, on trouve le programme *deformedShape* qui permet,
en grandes transformations, de visualiser les champs sorties par le programme dans un 
VTK pour lequel les voxels sont dans la configuration actuelle. Ce programme s'exécute simplement
en lançant la commande :

    ./deformedShape racine_in [racine_out]  
	
Le premier paramètre représente la racine du (des) fichier(s) d'entrée et
le deuxième paramètre est optionnel et permet de définir la racine du fichier de sortie.
Par défaut, cette racine est *sortie*.
Le programme va alors chercher le fichier *racine_in_def.vtk* ou, s'il ne le trouve pas,
les fichiers *racine_in_defi.vtk* (*i* allant de 1 à 9).
Ce champ représentant le gradient du déplacement permet alors de construire
le champ de déplacement et de calculer les coordonnées
de chaque voxel dans la configuration déformée. Le gradient du déplacement est alors
réécrit dans  cette configuration dans le fichier *racine_out.vtk*.
De plus, si les fichiers *racine_in_sig* ou *racine_in_pi*
(représentant respectivement les champs de contraintes de Cauchy et de 
Piola-Kirchhoff) sont présents, les variables correspondantes seront également présentes
dans le fichier *racine_out.vtk*.

5 - Pour les développeurs : validation et mise à jour des résultats
-------------------------------------------------------------------
Dans le dossier *validation*, le script *script_tests.sh* lance certains calculs
et vérifie que les résultats restent cohérents par rapport à ceux obtenus précédemment.
Le script *update_results.sh* met à jour les résultats de la base en prenant ceux
obtenus avec la nouvelle version du programme.

IMPORTANT
---------
Avant chaque push, le développeur devra lancer les tests pour vérifier que les resultats
obtenus sont justes et que les développements n'entraînent pas une regression du code.
Ces tests sont lancés avec un code compilé en utilisant intel15.
Il faudra ensuite mettre à jour les résultats des tests en lançant le script *update_results*.
Avant chaque commit, il est préférable de lancer le script clean pour éviter d'envoyer des chemins d'accès locaux.
