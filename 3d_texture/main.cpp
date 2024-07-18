/// Ce programme sert tout d'abord à reconstruire la cartographie des noyaux et spots, tester si tel spot appartient à tel noyau, présenter la cartographie de noyaux en fonction de nombre
/// de spots présents dans chaque noyau en donnant une taille différente ainsi qu'un intensité pixel différent pour chaque noyau.
/// Le nombre de noyaux et spots sont enregistrés dans nombre_cell.csv et nombre_spot.csv qui sont dans le dossier codeblocks, la cartographie (cellSfr_i.tif) y est aussi.
/// On peut classifier le nombre de spots sous différentes couches (10 couches dans ce programme), chaque couche présente une distance 0.1R par rapport au centre de sphéroide.
/// Le résultat est enregistré dans le fichier hist_spot.csv qui se trouve dans le dossier codeblocks.

/// 4 paramètres à modifier en fonction des expériences: userName (ligne 36), acquiDate (ligne 37), frame (ligne 38) et N3D(ligne 33, dimx, dimy, dimz, par défaut 1200 si on crop les images à 1200 avec macro crop_bleachcorrect.ijm)

#include "fonction.h"
#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <tiffio.h>
#include <math.h>

//#include <highgui.h>
//#include <cv.h>
//using namespace cv;

#ifndef DEF_STRUCT// Si la constante n'a pas été définie` le fichier n'a jamais été inclus
#define DEF_STRUCT

#endif // DEF_STRUCT

using namespace std;

int main()
{
    unsigned short dimx = 1200, dimy = 1200, dimz = 150;
    //unsigned short dimx = 1024, dimy = 1024, dimz = 150;
    Var3D N3D={dimx, dimy, dimz};
    string userName = "Corinne";
    string acquiDate = "170906";
    int frame = 50;
    lire_csv(userName, acquiDate, frame, N3D);
}

