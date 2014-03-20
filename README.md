AGGP
====

## Générateur de graphes biologiquement plausibles par algorithme génétique.

* Compte github à créer

Dans user settings, ssh keys : coller sa clef rsa, obtenue avec : 
Si plusieurs ordis, plusieurs clefs à mettre ;)

>ssh-keygen -t rsa -C "votremailidentifiant"
>cat ~/.ssh/id_rsa.pub

Me donner votre identifiant pour que je vous ajoute au projet. 



* Pour récupérer le code : première fois :

dans un terminal avec git d'installer, sinon sous Ubuntu.debian : apt-get install git

>$ git clone https://github.com/Fritzip/AGGP.git


ça crée le dossier avec les fichiers,

vous pouvez maintenant tout modifier normalement,

* Une fois le dossier créer avant de travailler, mettre à jour le dépot (pour prendre en compte les modifs des autres)

>$ git pull

(il faut etre dans le dossier)


* Après chaque changement IMPORTANT : faire un commit :

>$ git commit -a -m 'un message explicite !'

* A la fin de la session de travail : 

>$ git push -u origin master

cela met à jour le code sur le serveur, avant ce n'est que sur votre pc donc on en profite pas !

* pour ajouter un fichier :
 
>$ git add fichier

* pour voir l'état de git : (ce qui est modif ou pas, ce qui est inclu dans le dépot ou pas )

>$ git status


Un peu plus d'aide là, si vous voulez :

http://doc.ubuntu-fr.org

#### Merci à Antoine pour ce README.