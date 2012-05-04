#!/bin/sh
cd ..
cp update.pyw __main__.py
zip -r tools/install.pyw __main__.py src
rm __main__.py
