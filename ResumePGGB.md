# Analyse Méthodologique : Construction de Graphes de Pangenome (Garrison et al., 2024)

## Résumé Scientifique

Les graphes de pangenome constituent une représentation compacte et intégrée des génomes complets, permettant de capturer simultanément les homologies, les variations structurelles et les polymorphismes nucléotidiques entre génomes apparentés. Cette approche transcende les limitations inhérentes aux méthodes conventionnelles d'alignement de novo et de cartographie sur génome de référence unique, qui induisent des biais de représentation et conduisent à des pertes d'information génomique critique.

Face aux insuffisances des approches guidées par référence, susceptibles de générer des représentations incomplètes et instables de la variation génétique, Garrison et collaborateurs (2024) ont développé le PanGenome Graph Builder (PGGB), un pipeline modulaire implémenté en scripts shell et initialement validé sur le génome humain dans le cadre du Human Pangenome Reference Consortium (HPRC).

## Architecture Méthodologique

### Principe Fondamental : Approche Sans Référence

PGGB adopte une stratégie radicalement différente en traitant tous les génomes de manière équivalente, sans présupposition concernant les relations phylogénétiques, les groupes d'orthologie ou l'histoire évolutionnaire. Cette philosophie méthodologique s'articule autour d'alignements tous-contre-tous (all-to-all), où chaque séquence du pangenome peut potentiellement servir de référence pour décrire la variation.

**Complexité computationnelle :** Cette approche présente une complexité O(n²) pour les comparaisons, avec une croissance exponentielle de l'utilisation mémoire et un temps proportionnel au nombre de paires de génomes analysés.

### Phase I : Alignement par Segmentation (WFMASH)

**Innovation technique :** PGGB contourne les limitations computationnelles par l'utilisation de WFMASH, qui effectue la cartographie et l'alignement dans l'espace des segments de séquences plutôt qu'au niveau des paires de bases individuelles.

**Avantages :**
- Réduction drastique des coûts computationnels
- Préservation de la colinéarité à long terme
- Insensibilité aux similarités courtes répétitives (transposons, séquences satellites)

**Limitation :** Bien que WFMASH génère des alignements optimaux pour le fonctionnement de PGGB, le pipeline conserve la flexibilité d'accepter tout ensemble d'alignements définis par l'utilisateur au format PAF.

### Phase II : Induction du Graphe (SEQWISH)

**Principe :** SEQWISH convertit l'ensemble des génomes et alignements par paires en un graphe de variation équivalent, en fusionnant toutes les paires de bases appariées dans les alignements en un nœud unique dans le graphe résultant.

**Propriétés remarquables :**
- Compression des paires de bases transitivement appariées
- Récupération des relations d'homologie transitives absentes de l'ensemble d'alignements initial
- Permettant l'application de sparsification aléatoire pour réduire la complexité

**Stratégie de sparsification :** PGGB implémente une heuristique basée sur le modèle d'Erdős–Rényi pour déterminer un seuil sécurisé de sparsification. Le paramètre utilise un hash de chaque enregistrement de cartographie pour éliminer les mappings avec une probabilité Psparse ≫ Pconnected, garantissant la préservation des composantes connexes englobant les relations homologues critiques.

**Impact :** Cette approche évite les coûts O(N²) attendus lorsque Psparse=1, réduisant dramatiquement le temps d'exécution de l'alignement et de l'induction du graphe avec un effet négligeable sur la précision.

### Phase III : Normalisation du Graphe (SMOOTHXG)

**Problématique :** Bien que le graphe SEQWISH présente un modèle complet et sans perte des génomes d'entrée et de leurs homologies, il exhibe souvent des motifs locaux complexes problématiques pour les analyses en aval. Un enjeu majeur réside dans la non-normalisation mutuelle des alignements par paires, conduisant à des représentations différentielles des petites variations (indels) dans les séquences de faible complexité.

**Solution :** SMOOTHXG constitue une étape de post-traitement itératif spécifiquement conçue pour compresser et simplifier localement le graphe de variation du pangenome.

**Mécanisme :** Application d'un noyau de réalignement par Alignement d'Ordre Partiel (POA) aux régions extraites d'un embedding graphique unidimensionnel. Cet embedding ordonne les nœuds de sorte que leur distance dans l'ordre approxime au mieux leur distance dans les chemins génomiques du graphe.

**Processus itératif :**
1. Apprentissage de l'embedding
2. Obtention de segments partiellement chevauchants (blocs)
3. Application du POA aux blocs
4. "Laçage" des blocs réalignés dans un graphe de variation complet
5. Itération multiple (3 par défaut) pour limiter les effets de bord aux frontières des blocs

**Normalisation finale :** Application de GFAFFIX pour compresser les nœuds redondants et utilisation d'ODGI pour un tri final du graphe modifié.

## Validation et Outputs

### Contrôle Qualité Intégré

PGGB génère des outputs supportant l'interprétation immédiate, le contrôle qualité et les applications en aval :
- Statistiques basiques du graphe (taille, nombre de nœuds, contenu en bases)
- Visualisations 1D et 2D via ODGI
- Rapports diagnostiques optionnels via MultiQC

### Appel de Variants

**Innovation :** PGGB fonctionne comme un appeleur de variants multi-échantillons pour assemblages de génomes complets, générant une description phasée des haplotypes embarqués au format VCF.

**Défi technique :** Les variants appelés directement du graphe peuvent inclure des sites génétiques grands et imbriqués, créant des incompatibilités avec de nombreuses applications.

**Solution :** Décomposition de la variation complexe imbriquée en une représentation minimale relative à la référence utilisant BiWFA, assurant la compatibilité avec les analyses basées sur les petites variations.

## Validation Croisée et Performance

### Validation par MUMMER4

La validation croisée avec MUMMER4, standard actuel pour l'alignement génomique par paires, démontre des F-scores >92% dans tous les contextes testés, indiquant une performance équivalente aux standards existants.

**Divergences observées :** Les différences entre PGGB et MUMMER4 résultent de représentations d'alignement distinctes dans les régions présentant des SNVs très proches.

**Avantage comparatif :** Contrairement à MUMMER4 qui ne fournit que des comparaisons par paires avec une référence cible, PGGB génère une comparaison complète tous-contre-tous permettant de nouvelles modalités d'analyse bioinformatique.

### Applications Phylogénétiques

**Preuve de concept :** Construction d'arbres phylogénétiques directement à partir de distances mesurées dans un graphe de variation de pangenome de 16 assemblages complets du chromosome 6 de la famille des grands singes, concordant avec les phylogénies établies du clade Hominoidea.

**Validation comparative :** Application réussie chez la levure, où la comparaison d'arbres phylogénétiques construits à partir des SNVs de MUMMER contre la référence et des nœuds du graphe PGGB montre que ce dernier infère correctement la phylogénie dans les cas où les deux haplotypes du même génome sont rapprochés.

## Limitations et Perspectives

### Limitations Identifiées

1. **Complexité computationnelle :** Bien que réduite par WFMASH, l'approche reste quadratique
2. **Partitionnement requis :** Les génomes volumineux nécessitent un partitionnement pour maximiser le parallélisme
3. **Représentation des variants complexes :** Potentielles incompatibilités avec certaines applications en aval
4. **Absence de vérité terrain :** Difficultés d'évaluation qualitative sur données réelles

### Forces Méthodologiques

1. **Absence de biais de référence :** Traitement équitable de tous les génomes
2. **Scalabilité :** Performance sur des centaines de petits génomes en quelques heures
3. **Modularité :** Architecture flexible et extensible
4. **Exhaustivité :** Capture de toutes les classes de variation simultanément

## Conclusion Méthodologique

PGGB représente une approche conceptuellement novatrice et méthodologiquement robuste pour comprendre les relations séquentielles entre génomes complets multiples, tant dans les contextes pangénomiques que de génomique comparative. L'architecture modulaire et la philosophie sans référence ouvrent la voie à des méthodes génétiques populationnelles et évolutives diversifiées, capables de considérer simultanément toutes les classes de variation séquentielle.

Cette approche révolutionnaire permet le développement d'une compréhension exhaustive des liens entre variation séquentielle, phénotype et évolution dans une ère où l'assemblage complet des génomes devient routinier.

---

**Note méthodologique :** L'implémentation Nextflow scalable en cluster est actuellement en développement, renforçant l'applicabilité du pipeline aux projets à grande échelle.