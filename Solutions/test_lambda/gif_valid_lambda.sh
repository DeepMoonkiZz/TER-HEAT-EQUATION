#!/bin/bash

gnuplot gif_solution_valid_lambda.txt

# Définition des paramètres
n_files=101
t0=0.000
tfinal=0.100
dt=$(awk "BEGIN {print ($tfinal - $t0) / ($n_files - 1)}")

# Chemin d'accès au répertoire contenant les fichiers PNG
path="$PWD/validation_lambda"

# Création des répertoires pour les groupes d'images
mkdir -p group1 group2 group3

# Déplacer les images dans les répertoires appropriés
for i in $(seq 0 $((n_files - 1))); do
    t=$(awk "BEGIN {printf \"%.3f\", $t0 + $i * $dt}" | tr ',' '.')
    file_name="valid_lambda_$t.png"
    
    if [[ $i -lt 30 ]]; then
        mv "$path/$file_name" group1/
    elif [[ $i -lt 60 ]]; then
        mv "$path/$file_name" group2/
    else
        mv "$path/$file_name" group3/
    fi
done

# Générer les GIF pour chaque groupe avec une taille limitée à 800 pixels de largeur
convert -delay 10 -loop 0 group1/valid_lambda_*.png group1.gif
convert -delay 10 -loop 0 group2/valid_lambda_*.png group2.gif
convert -delay 10 -loop 0 group3/valid_lambda_*.png group3.gif

# Combinaison des GIF en un seul
gifsicle --merge group1.gif group2.gif group3.gif > /home/segal/Documents/MatMeca/TER-2A/Solutions/gif/valid_lambda.gif

# Supprimer les répertoires temporaires
rm -r group1 group2 group3
rm -f group1.gif group2.gif group3.gif
rm -f validation_lambda/*.png
