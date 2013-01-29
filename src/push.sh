git commit --interactive
while git status|grep "modified:" > /dev/null; do
    git commit --interactive
done && git push
git checkout master
git merge experimental
git push
git checkout experimental
rsync -ar --no-owner --no-g --modify-window=5 --delete ./ ~/documentos/carlos/Dropbox/golsoft/src

cd ~/documentos/carlos/Dropbox/golsoft/src

#for file in *.py *.c *.pyw; do
#    todos "$file"
#done
