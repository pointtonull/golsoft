#!/bin/sh
echo "Instalando nuevas dependencias"
sudo aptitude install -y python-gdal mayavi2 python-numpy python-matplotlib python-pygame python-opencv python-imaging
cd ~/Dropbox/golsoft/src
python run_pea_traits.pyw
