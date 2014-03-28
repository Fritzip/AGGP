AGGP
====

#### Générateur de graphes biologiquement plausibles par algorithme génétique.

* Installer Git

`apt-get install git`

* Créer un compte GitHub

**Dans `Acount settings -> SSH keys` coller sa clef RSA, obtenue avec : 

```ssh-keygen -t rsa -C \<votremailidentifiant\> 

cat ~/.ssh/id_rsa.pub```

Si plusieurs ordis, plusieurs clefs à mettre ;)

**Me donner votre identifiant pour que je vous ajoute au projet. 



* Pour récupérer le code pour la première fois :

`git clone https://github.com/Fritzip/AGGP.git`

ça crée le dossier AGGP avec les fichiers, vous pouvez maintenant tout modifier normalement,

* Une fois le dossier créer avant *chaque session de travail*, mettre à jour le dépot (pour prendre en compte les modifs des autres)

`git pull`

(il faut etre dans le dossier)


* Après chaque changement IMPORTANT, faire un commit :

`git commit -a -m 'un message explicite !'`

* A la fin de la session de travail (et si tout fonctionne bien !), envoyer sur le dépot:

`git push -u origin master`

cela met à jour le code sur le serveur, avant ce n'est que sur votre pc donc on en profite pas !

* Pour ajouter un fichier:
 
`git add \<fichier\>`

* Pour voir l'état de git (ce qui est modif ou pas, ce qui est inclu dans le dépot ou pas ) :

`git status`


Un peu [plus d'aide là](http://doc.ubuntu-fr.org), si vous voulez :


*Merci à Antoine pour ce README.*
