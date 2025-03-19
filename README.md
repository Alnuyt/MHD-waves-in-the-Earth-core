# Modèle simplifié d’ondes magnétohydrodynamiques dans le noyau terrestre
Dans ce mémoire, les équations linéaires de l’hydrodynamique et de la magnéto-
hydrodynamique sont développées pour un fluide conducteur en rotation uniforme
imprégné d’un champ magnétique constant. Une décomposition toroïdale-poloïdale
est apppliquée aux champs solénoïdaux, ce qui nous permet, d’obtenir un système
résoluble numériquement. La résolution numérique de ces équations via un solveur
d’ondes inertielles, développé dans le cadre de ce mémoire, permet de calculer tous
les modes d’oscillations de ces ondes pour le cas hydrodynamique. La base d’un
tel solveur a été construite durant cette étude, des améliorations mèneront vers un
véritable calculateur d’ondes magnéto-Coriolis applicable à la géophysique et l’as-
trophysique.
## Auteur
Alexandre Nuyt - Université Catholique de Louvain, Faculté des sciences, Ecole de physique
## Objectifs
Ce mémoire étant exploratoire, trois objectifs majeurs sont visés :
- Étudier les fondements de la dynamique des fluides en rotation et de la ma-
gnétohydrodynamique, dans un contexte géophysique.
- Développer les équations différentielles des ondes inertielles pour l’hydrody-
namique et la MHD en décomposition Toroïdale-Poloïdale.
- Construire la base d'un outil numérique : solveur Hydro + MHD des ondes inertielles,
pour tout mode m et k, en géométrie cylindrique.
## Description des codes
Lors de l'étude, 3 codes principaux ont été écrits :
- [Analytical_Inertial_Waves.ipynb](Analytical_Inertial_Waves.ipynb) : Résoud numériquement l'équation transcendante (1.56) qui donne les solutions analytiques des ondes inertielles se propageant dans fluide non-conducteur en rotation pour une géométrie cylindrique. 
- [Helmoltz_cylindrical.py](Helmoltz_cylindrical.py) : Résoud l'équation de Helmoltz cylindrique, par la méthode des diffénreces finies, afin de de s'assurer du fonctionnement de la méthode numérique utilisée pour réaliser les solveurs Hydro et MHD des ondes inertielles.
- [Toroidal_Poloidal_Inertial_waves.py](Toroidal_Poloidal_Inertial_waves.py): Est la base d'un solveur Hydro + MHD des ondes inertielles, pour tout mode m et k, en géométrie cylindrique.
## Conclusion
Ce travail explore la résolution des ondes magnéto-Coriolis en géométrie cylindrique en s'appuyant sur des bases théoriques en hydrodynamique et magnétohydrodynamique. Un algorithme numérique a été développé pour traiter ce problème, bien que les résultats actuels nécessitent encore des améliorations. En optimisant la discrétisation, cet outil pourra devenir un véritable atout pour l'étude des phénomènes géophysiques et astrophysiques.
