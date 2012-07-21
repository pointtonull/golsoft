rsync -ar --no-owner --no-g --modify-window=5 ~/documentos/carlos/Dropbox/golsoft/src/ ./

for file in *.py *.pyw *.sh *.c
    fromdos "$file"
