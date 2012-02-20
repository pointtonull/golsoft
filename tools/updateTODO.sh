#!/bin/sh
NOTEBOOK="~/documentos/carlos/apuntes/zim"
TODO="Home:Seminario:TODO"
DONE="Home:Seminario:TODO:DONE"
FULL="Home:Seminario:TODO:FULL"

zim --export --template Presentation --output=docs "$NOTEBOOK" "$TODO"
zim --export --template Presentation --output=docs "$NOTEBOOK" "$DONE"
zim --export --template Presentation --output=docs "$NOTEBOOK" "$FULL"

awk '{gsub("file:///usr/share/zim/", "");
gsub("./TODO/", "./TODO.")}/./' docs/TODO.html > docs/TODO.tmp
mv docs/TODO.tmp docs/TODO.html

awk '{gsub("file:///usr/share/zim/", "")}/./' docs/DONE.html > docs/DONE.tmp
mv docs/DONE.tmp docs/TODO.DONE.html
rm docs/DONE.html

awk '{gsub("file:///usr/share/zim/", "")}/./' docs/FULL.html > docs/FULL.tmp
mv docs/FULL.tmp docs/TODO.FULL.html
rm docs/FULL.html

ln -f ~/documentos/carlos/seminario/documento/Maestro.pdf docs/seminario\ Carlos.pdf

git add docs/*.html docs/*.pdf
git commit -m "updated docs"
git push origin master
