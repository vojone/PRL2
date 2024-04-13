#preklad zdrojoveho souboru
mpic++ --prefix /usr/local/share/OpenMPI -o life life.cpp

#spusteni programu
mpirun --prefix /usr/local/share/OpenMPI  -np $1 life			

#uklid
rm -f life
