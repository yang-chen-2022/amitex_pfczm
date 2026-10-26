
METHODE DE TRAVAIL POUR INTRODUIRE DES DEVELOPPEMENTS :
=====================================================

1. Cibler le développement a realiser et son type (fonctionnalité, debug), et 
ouvrir une **Issue** sur le Gitlab.

2. Créer une **Merge Request** sur la page de l'**Issue**. Celà crée automatiquement
une branche avec pour nom *NumeroIssue-Titre-de-l'-Issue*.

3. Récupérer la branche dédiée et implémenter la solution et les tests dessus.
```shell
git checkout --track remotes/origin/NumeroIssue-Titre-de-l'-Issue
```

4. Une fois les développements terminés, réaliser un premier *merge* de la 
branche *master* dans la branche *NumeroIssue-Titre-de-l'-Issue*.

5. Prévenir un **Reviewer** pour le développement et lui assigner l'**Issue** sur
Gitlab (en général, Lionel Gélébart). 

6. Le **Reviewer** passe les développements en revue et valide ou non le merge.

7. Suppression de la branche *NumeroIssue-Titre-de-l'-Issue*.
