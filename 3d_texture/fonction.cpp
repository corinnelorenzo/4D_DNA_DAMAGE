#include "fonction.h"

std::istream& operator>>(std::istream& str, CSVRow& data)
{
    data.readNextRow(str);
    return str;
}

void Read_bin(std::string chemin, unsigned char resultat[], short precision, int NbPix)
{
    long lTaille;
    int nb_elmnt_lu;
    unsigned short dimData=precision/8;//taille en octet d'un element.

    FILE* pFichier = NULL;

    pFichier = fopen(chemin.c_str(), "r");  //ouverture de ce fichier en écriture binaire
    std::cout<<pFichier<<std::endl;

    if(pFichier==NULL)
    {
        fputs("Impossible d'ouvrir le fichier\n",stderr);
        exit (1);// obtenir la longueur du fichier, comparer avec donnée entrée.
    }
    else
    {
        fseek(pFichier,0,SEEK_END);//trouver la fin de fichier

        lTaille = ftell (pFichier);//retourne la position courante (en octet) du curseur de fichier : ici, position de la fin du fichier
        printf("taille trouvée en octet par ftell %li, taille estimée : %i\n",lTaille,NbPix*dimData);//
        rewind(pFichier);

        if(NbPix*dimData!=lTaille)
        std::cout<<"Taille du fichier incompatible avec les dimensions\n"<<std::endl;

        nb_elmnt_lu = fread (resultat,1,lTaille,pFichier);//lecture

        if(nb_elmnt_lu!=lTaille){
            std::cout<<"Problème lors de la lecture du fichier"<<std::endl;

            std::cout<<"Nombre d'éléments lus="<<nb_elmnt_lu<<std::endl;
        }

        fclose(pFichier);
    }
}

void Write_bin_uchar(unsigned char *var_sav, int NbPix2D, std::string chemin, char options[])
{
        FILE *fichier_ID;
        fichier_ID= fopen(chemin.c_str(), options);
        if(fichier_ID==0)
            std::cout<<"Erreur d'ouverture du fichier "<<chemin<<std::endl;

        ///8 bit
        for(int cpt=0; cpt<NbPix2D; cpt++) {
            unsigned char tampon=var_sav[cpt];
            fwrite(&tampon,sizeof(tampon),1,fichier_ID);
        }
        fclose(fichier_ID);
}

void Write_bin_ushort(unsigned short *var_sav, int NbPix2D, std::string chemin, char options[])
{
        FILE *fichier_ID;
        fichier_ID= fopen(chemin.c_str(), options);
        if(fichier_ID==0)
            std::cout<<"Erreur d'ouverture du fichier "<<chemin<<std::endl;

        ///16 bit non signé
        for(int cpt=0; cpt<NbPix2D; cpt++) {
            unsigned short tampon=var_sav[cpt];
            fwrite(&tampon,sizeof(tampon),1,fichier_ID);
        }
        fclose(fichier_ID);
}

void WRITE_tiff(unsigned short *var_sav, std::string chemin, int dim)
{
    const char *chem = chemin.c_str();
    TIFF *tif= TIFFOpen(chem, "a");
    TIFFSetField (tif, TIFFTAG_IMAGEWIDTH, dim);
    TIFFSetField (tif, TIFFTAG_IMAGELENGTH, dim);
    TIFFSetField (tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField (tif, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField (tif, TIFFTAG_ROWSPERSTRIP, 1);
    TIFFSetField (tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
    TIFFSetField (tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    TIFFSetField (tif, TIFFTAG_BITSPERSAMPLE, 16);
    TIFFSetField (tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);
    TIFFSetField (tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);

    tsize_t strip_size = TIFFStripSize (tif);
    tstrip_t strips_num = TIFFNumberOfStrips (tif);

    unsigned short* strip_buf=(unsigned short*)_TIFFmalloc(strip_size);
    for (int s=0; s<strips_num; s++)
    {
        for (int col=0; col<dim; col++)
        {
            int cpt=col+dim*s;
            strip_buf[col]=var_sav[cpt];
        }
        TIFFWriteEncodedStrip (tif, s, strip_buf, strip_size);
    }
    _TIFFfree(strip_buf);
    TIFFWriteDirectory(tif);
    TIFFClose(tif);
}

void WRITE_tiff3D(unsigned short *var_sav, std::string chemin, Var3D dim)
{
    unsigned short *buffer=new unsigned short[dim.x*dim.y];
    for(int z=0; z < dim.z; z++)
    {
        for(int y=0; y < dim.y; y++)
        {
            for(int x=0; x < dim.x; x++)
            {
                int cpt1= x + y * dim.x + z * dim.x * dim.y;
                int cpt2= x + y * dim.x;
                buffer[cpt2]=var_sav[cpt1];
            }
        }
        WRITE_tiff(buffer, chemin, dim.x);
	    }
    delete[] buffer;
}

int maximum( int a, int b, int c )
{
   int max = ( a < b ) ? b : a;
   return ( ( max < c ) ? c : max );
}

//void genere_SFR(float *objet, Var3D cObj, Var3D rObj, int R_max, int mInt, Var3D dim, int i)
//{
//    double Rcarre = R_max*R_max;
//    //double Rcarre = (rObj.x*rObj.x + rObj.y*rObj.y + rObj.z*rObj.z)/3;
//    for(int z = cObj.z-rObj.z; z < cObj.z+rObj.z; z++){
//        if ((z+rObj.z)<=dim.z && z>0) /// segmentation fault si z+R dépasse dim.z (150)
//        {
//        int altitude = dim.x*dim.y*z;
//        for(int y = cObj.y-rObj.y; y < cObj.y+rObj.y; y++){
//            int NbPixY = dim.x*y;
//            for(int x = cObj.x-rObj.x; x<cObj.x+rObj.x; x++){
//                double d2 = pow((x-cObj.x),2)+pow((y-cObj.y),2)+pow((z-cObj.z),2);
//                if(d2<Rcarre){
//                    int cpt3D = altitude+NbPixY+x;
//                    //objet[cpt3D] = mInt;   /// valeur = mean intensity
//                    objet[cpt3D] = i;   /// valeur = numéro de bbox
//                }
//                else{
//                    int cpt3D = altitude+NbPixY+x;
//                    objet[cpt3D] = 0;
//                }
//
//            }
//        }
//        }
//    }
//}

void genere_rayonfix_cell(unsigned short *objet, Var3D cObj, Var3D rBox, int Rfixe, Var3D dim, int i)
{
    double Rcarre = Rfixe*Rfixe;
    //double Rcarre = (rObj.x*rObj.x + rObj.y*rObj.y + rObj.z*rObj.z)/3;

    int bDim = rBox.x * 2 * rBox.y * 2 * rBox.z * 2;
    if ((bDim * 0.32 * 0.32 < 6000) && (bDim * 0.32 * 0.32 > 100)){  /// Cell
    for(int z = cObj.z-Rfixe; z < cObj.z+Rfixe; z++){
        if ((z+Rfixe)<=dim.z && z>0) /// segmentation fault si z+R dépasse dim.z (150)
        {
        int altitude = dim.x*dim.y*z;
        for(int y = cObj.y-Rfixe; y < cObj.y+Rfixe; y++){
            int NbPixY = dim.x*y;
            for(int x = cObj.x-Rfixe; x<cObj.x+Rfixe; x++){
                double d2 = pow((x-cObj.x),2)+pow((y-cObj.y),2)+pow((z-cObj.z),2);
                if(d2<Rcarre){
                    int cpt3D = altitude+NbPixY+x;
                    objet[cpt3D] = i;   ///
                }
                else{
                    int cpt3D = altitude+NbPixY+x;
                    objet[cpt3D] = 0;
                }
            }
        }
        }
    }
    }
}

void genere_SFR_cell(unsigned short *objet, Var3D cObj, Var3D rBox, int Rbox_max, int mInt, Var3D dim, int i)
{
    double Rcarre = Rbox_max*Rbox_max;
    //double Rcarre = (rObj.x*rObj.x + rObj.y*rObj.y + rObj.z*rObj.z)/3;

    int bDim = rBox.x * 2 * rBox.y * 2 * rBox.z * 2;
    if ((bDim * 0.32 * 0.32 < 6000) && (bDim * 0.32 * 0.32 > 100)){  /// Cell
    for(int z = cObj.z-rBox.z; z < cObj.z+rBox.z; z++){
        if ((z+rBox.z)<=dim.z && z>0) /// segmentation fault si z+R dépasse dim.z (150)
        {
        int altitude = dim.x*dim.y*z;
        for(int y = cObj.y-rBox.y; y < cObj.y+rBox.y; y++){
            int NbPixY = dim.x*y;
            for(int x = cObj.x-rBox.x; x<cObj.x+rBox.x; x++){
                double d2 = pow((x-cObj.x),2)+pow((y-cObj.y),2)+pow((z-cObj.z),2);
                if(d2<Rcarre){
                    int cpt3D = altitude+NbPixY+x;
                    //objet[cpt3D] = mInt;   /// valeur = mean intensity
                    objet[cpt3D] = i;   /// valeur = numéro de bbox
                }
                else{
                    int cpt3D = altitude+NbPixY+x;
                    objet[cpt3D] = 0;
                }
            }
        }
        }
    }
    }
}

void genere_BBox_cell(unsigned short *objet, Var3D bMin, Var3D bMax, Var3D cObj, int mInt, int sObj, Var3D dim, int i)
{
    int bDim = (bMax.z - bMin.z) * (bMax.y - bMin.y) * (bMax.x - bMin.x);
    //if ((bDim * 0.32 * 0.32 < 6000) && (bDim * 0.32 * 0.32 > 100)){  /// Cell
    if ((bDim * 0.32 * 0.32 < 10000) && (bDim * 0.32 * 0.32 > 100)){  /// Cell
        for(int z = bMin.z; z <= bMax.z; z++){
            int altitude = dim.x*dim.y*z;
            for(int y = bMin.y; y <= bMax.y; y++){
                int NbPixY = dim.x*y;
                for(int x = bMin.x; x <= bMax.x; x++){
                    int cpt3D = altitude+NbPixY+x;
                    //if(abs((z-cObj.z)*(y-cObj.y)*(x-cObj.x)) < sObj)
                    //{
                    objet[cpt3D] = i;   /// valeur = numéro de bbox
                    //objet[cpt3D] = mInt;   /// valeur = mean intensity
                    }
                    //}
            }
        }
    }
}

void genere_SFR_spot(unsigned short *objet, Var3D cObj, Var3D rBox, int Rbox_max, int mInt, Var3D dim, int i)
{
    double Rcarre = Rbox_max*Rbox_max;
    //double Rcarre = (rObj.x*rObj.x + rObj.y*rObj.y + rObj.z*rObj.z)/3;

    int bDim = rBox.x * 2 * rBox.y * 2 * rBox.z * 2;
    if ((bDim * 0.32 * 0.32 < 100) && (bDim * 0.32 * 0.32 > 0)){  /// spot

    for(int z = cObj.z-rBox.z; z < cObj.z+rBox.z; z++){
        if ((z+rBox.z)<=dim.z && z>0) /// segmentation fault si z+R dépasse dim.z (150)
        {
        int altitude = dim.x*dim.y*z;
        for(int y = cObj.y-rBox.y; y < cObj.y+rBox.y; y++){
            int NbPixY = dim.x*y;
            for(int x = cObj.x-rBox.x; x<cObj.x+rBox.x; x++){
                double d2 = pow((x-cObj.x),2)+pow((y-cObj.y),2)+pow((z-cObj.z),2);
                if(d2<Rcarre){
                    int cpt3D = altitude+NbPixY+x;
                    //objet[cpt3D] = mInt;   /// valeur = mean intensity
                    objet[cpt3D] = i;   /// valeur = numéro de bbox
                }
                else{
                    int cpt3D = altitude+NbPixY+x;
                    objet[cpt3D] = 0;
                }

            }
        }
        }
    }
    }
}

void genere_BBox_spot(unsigned short *objet, Var3D bMin, Var3D bMax, Var3D cObj, int mInt, int sObj, Var3D dim, int i)
{
    int bDim = (bMax.z - bMin.z) * (bMax.y - bMin.y) * (bMax.x - bMin.x);
    if ((bDim * 0.32 * 0.32 < 100) && (bDim * 0.32 * 0.32 > 0)){  /// spot
        for(int z = bMin.z; z <= bMax.z; z++){
            int altitude = dim.x*dim.y*z;
            for(int y = bMin.y; y <= bMax.y; y++){
                int NbPixY = dim.x*y;
                for(int x = bMin.x; x <= bMax.x; x++){
                    int cpt3D = altitude+NbPixY+x;
                    //if(abs((z-cObj.z)*(y-cObj.y)*(x-cObj.x)) < sObj)
                    //{
                    objet[cpt3D] = i;   /// valeur = numéro de bbox
                    //objet[cpt3D] = mInt;   /// valeur = mean intensity
                    }
                    //}
            }
        }
    }
}

void lire_csv(std::string userName, std::string acquiDate, int frame, Var3D N3D)
{
    /// lire fichier.csv (généré avec Ilastik)
    ///cell
    const int N = 100000;
    std::vector<float> sizeObject_cell(N);
    std::vector<float> bbMin_x_cell(N);
    std::vector<float> bbMin_y_cell(N);
    std::vector<float> bbMin_z_cell(N);
    std::vector<float> bbMax_x_cell(N);
    std::vector<float> bbMax_y_cell(N);
    std::vector<float> bbMax_z_cell(N);
    std::vector<float> meanIntensity_cell(N);
    std::vector<float> centerObject_x_cell(N);
    std::vector<float> centerObject_y_cell(N);
    std::vector<float> centerObject_z_cell(N);
    std::vector<float> rayon_x_cell(N);
    std::vector<float> rayon_y_cell(N);
    std::vector<float> rayon_z_cell(N);

    ///spot
    std::vector<float> sizeObject_spot(N);
    std::vector<float> bbMin_x_spot(N);
    std::vector<float> bbMin_y_spot(N);
    std::vector<float> bbMin_z_spot(N);
    std::vector<float> bbMax_x_spot(N);
    std::vector<float> bbMax_y_spot(N);
    std::vector<float> bbMax_z_spot(N);
    std::vector<float> meanIntensity_spot(N);
    std::vector<float> centerObject_x_spot(N);
    std::vector<float> centerObject_y_spot(N);
    std::vector<float> centerObject_z_spot(N);
    std::vector<float> rayon_x_spot(N);
    std::vector<float> rayon_y_spot(N);
    std::vector<float> rayon_z_spot(N);

    std::ofstream hist_spot;
    std::string file_hist_spot = "E:\\" + userName + "\\" + acquiDate + "\\codeblocks\\hist_spot.csv";
    hist_spot.open(file_hist_spot.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);

    std::ofstream nombre_cell;
    std::string file_nombre_cell = "E:\\" + userName + "\\" + acquiDate + "\\codeblocks\\nombre_cell.csv";
    nombre_cell.open(file_nombre_cell.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);

    std::ofstream nombre_spot;
    std::string file_nombre_spot = "E:\\" + userName + "\\" + acquiDate + "\\codeblocks\\nombre_spot.csv";
    nombre_spot.open(file_nombre_spot.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);


    for(int t=1; t<=frame; t++)
    {
        unsigned short * cellSfr = new unsigned short [N3D.x*N3D.y*N3D.z];
        unsigned short * distSfr = new unsigned short [N3D.x*N3D.y*N3D.z];

        ///cell
        std::ostringstream oss;
        oss << t;
        printf("Temps = T%i\n", t);

        std::string file_cell = "E:\\" + userName + "\\" + acquiDate + "\\cell\\csv\\" + acquiDate +  "_N_" + oss.str() + ".csv"; ///N: noyau
        std::ifstream ifs_file_cell(file_cell.c_str());

        CSVRow  row_cell;
        int lineNumber_cell = 0;
        while(ifs_file_cell >> row_cell)
        {
            sizeObject_cell[lineNumber_cell] = atof(row_cell[7].c_str());

            bbMin_x_cell[lineNumber_cell] = atof(row_cell[8].c_str());
            bbMin_y_cell[lineNumber_cell] = atof(row_cell[9].c_str());
            bbMin_z_cell[lineNumber_cell] = atof(row_cell[10].c_str());

            centerObject_x_cell[lineNumber_cell] = atof(row_cell[11].c_str());
            centerObject_y_cell[lineNumber_cell] = atof(row_cell[12].c_str());
            centerObject_z_cell[lineNumber_cell] = atof(row_cell[13].c_str());

            bbMax_x_cell[lineNumber_cell] = atof(row_cell[14].c_str());
            bbMax_y_cell[lineNumber_cell] = atof(row_cell[15].c_str());
            bbMax_z_cell[lineNumber_cell] = atof(row_cell[16].c_str());

            rayon_x_cell[lineNumber_cell] = atof(row_cell[31].c_str());
            rayon_y_cell[lineNumber_cell] = atof(row_cell[32].c_str());
            rayon_z_cell[lineNumber_cell] = atof(row_cell[33].c_str());

            meanIntensity_cell[lineNumber_cell] = atof(row_cell[113].c_str());

            lineNumber_cell ++;
        }
        std::cout << "lineNumber_cell = " << lineNumber_cell << std::endl;
        nombre_cell << lineNumber_cell << std::endl;

//    for(int i = 0; i < lineNumber_cell; i ++)
//    {
//        printf("BBox_min_cell: x = %f, y = %f, z = %f, BBox_max_cell: x = %f, y = %f, z = %f\n", bbMin_x_cell[i], bbMin_y_cell[i], bbMin_z_cell[i], bbMax_x_cell[i], bbMax_y_cell[i], bbMax_z_cell[i]);
//        printf("Mean intensity_cell: %f, CenterObject_cell: x = %f, y = %f, z = %f \n", meanIntensity_cell[i], centerObject_x_cell[i], centerObject_y_cell[i], centerObject_z_cell[i]);
//    }

        ///spot
        std::string file_spot = "E:\\" + userName + "\\" + acquiDate + "\\spot\\csv\\" + acquiDate +  "_S_" + oss.str() + ".csv"; ///S: spot
        std::ifstream ifs_file_spot(file_spot.c_str());

        CSVRow  row_spot;
        int lineNumber_spot = 0;
        while(ifs_file_spot >> row_spot)
        {
            sizeObject_spot[lineNumber_spot] = atof(row_spot[7].c_str());

            bbMin_x_spot[lineNumber_spot] = atof(row_spot[8].c_str());
            bbMin_y_spot[lineNumber_spot] = atof(row_spot[9].c_str());
            bbMin_z_spot[lineNumber_spot] = atof(row_spot[10].c_str());

            centerObject_x_spot[lineNumber_spot] = atof(row_spot[11].c_str());
            centerObject_y_spot[lineNumber_spot] = atof(row_spot[12].c_str());
            centerObject_z_spot[lineNumber_spot] = atof(row_spot[13].c_str());

            bbMax_x_spot[lineNumber_spot] = atof(row_spot[14].c_str());
            bbMax_y_spot[lineNumber_spot] = atof(row_spot[15].c_str());
            bbMax_z_spot[lineNumber_spot] = atof(row_spot[16].c_str());

            rayon_x_spot[lineNumber_spot] = atof(row_spot[31].c_str());
            rayon_y_spot[lineNumber_spot] = atof(row_spot[32].c_str());
            rayon_z_spot[lineNumber_spot] = atof(row_spot[33].c_str());

            meanIntensity_spot[lineNumber_spot] = atof(row_spot[113].c_str());

            lineNumber_spot ++;
        }
        std::cout << "lineNumber_spot = " << lineNumber_spot << std::endl;
        nombre_spot << lineNumber_spot << std::endl;
    /*
    for(int i = 0; i < lineNumber_spot; i ++)
    {
        printf("BBox_min_spot: x = %f, y = %f, z = %f, BBox_max_spot: x = %f, y = %f, z = %f\n", bbMin_x_spot[i], bbMin_y_spot[i], bbMin_z_spot[i], bbMax_x_spot[i], bbMax_y_spot[i], bbMax_z_spot[i]);
        printf("Mean intensity_spot: %f, CenterObject_spot: x = %f, y = %f, z = %f \n", meanIntensity_spot[i], centerObject_x_spot[i], centerObject_y_spot[i], centerObject_z_spot[i]);
    }*/

/// début comment, si on veut juste savoir le nombre de noyaux et spots

    int cpt;
    //ofstream myfile ("spots.txt");

    int Cx_max=0, Cy_max=0, Cz_max=0;
    int Cx_min=N3D.x, Cy_min=N3D.y, Cz_min=N3D.z;

    for(int i = 1; i <= lineNumber_cell; i ++){

        int Cx_cell = round(centerObject_x_cell[i]);
        int Cy_cell = round(centerObject_y_cell[i]);
        int Cz_cell = round(centerObject_z_cell[i]);
        int Rx_cell = round(rayon_x_cell[i]);
        int Ry_cell = round(rayon_y_cell[i]);
        int Rz_cell = round(rayon_z_cell[i]);
        int bMinx_cell = round(bbMin_x_cell[i]);
        int bMiny_cell = round(bbMin_y_cell[i]);
        int bMinz_cell = round(bbMin_z_cell[i]);
        int bMaxx_cell = round(bbMax_x_cell[i]);
        int bMaxy_cell = round(bbMax_y_cell[i]);
        int bMaxz_cell = round(bbMax_z_cell[i]);

        Var3D cObj_cell = {Cx_cell, Cy_cell, Cz_cell};
        Var3D rObj_cell = {Rx_cell, Ry_cell, Rz_cell};
        int mInt_cell = round(meanIntensity_cell[i]);
        int sObj_cell = round(sizeObject_cell[i]);

        int rbox_x_cell = (bMaxx_cell -bMinx_cell)/2;
        int rbox_y_cell = (bMaxy_cell -bMiny_cell)/2;
        int rbox_z_cell = (bMaxz_cell -bMinz_cell)/2;
        int R_max_cell = maximum(Rx_cell, Ry_cell, Rz_cell);
        int Rbox_max_cell = maximum(rbox_x_cell, rbox_y_cell, rbox_z_cell);

        Var3D bMin_cell = {bMinx_cell, bMiny_cell, bMinz_cell};
        Var3D bMax_cell = {bMaxx_cell, bMaxy_cell, bMaxz_cell};
        Var3D rBox_cell = {rbox_x_cell, rbox_y_cell, rbox_z_cell};

        int bDim = rbox_x_cell * 2 * rbox_y_cell * 2 * rbox_z_cell * 2;
        /// calcul la position maximum et minimum des cellules afin de calculer le centre de sphéroide
        if (Cx_cell > Cx_max)
        {
            Cx_max = Cx_cell;
        }
        if (Cx_cell < Cx_min)
        {
           Cx_min = Cx_cell;
        }
        if (Cy_cell > Cy_max)
        {
            Cy_max = Cy_cell;
        }
        if (Cy_cell < Cy_min)
        {
           Cy_min = Cy_cell;
        }
        if (Cz_cell > Cz_max)
        {
            Cz_max = Cz_cell;
        }
        if (Cz_cell < Cz_min)
        {
           Cz_min = Cz_cell;
        }

        ///générer la carte des cellules en fonction de nombre de spots présents dans la cellule
        cpt = 0;
        for(int j = 0; j <= lineNumber_spot; j ++){

            int Cx_spot = round(centerObject_x_spot[j]);
            int Cy_spot = round(centerObject_y_spot[j]);
            int Cz_spot = round(centerObject_z_spot[j]);
            int Rx_spot = round(rayon_x_spot[j]);
            int Ry_spot = round(rayon_y_spot[j]);
            int Rz_spot = round(rayon_z_spot[j]);
            int bMinx_spot = round(bbMin_x_spot[j]);
            int bMiny_spot = round(bbMin_y_spot[j]);
            int bMinz_spot = round(bbMin_z_spot[j]);
            int bMaxx_spot = round(bbMax_x_spot[j]);
            int bMaxy_spot = round(bbMax_y_spot[j]);
            int bMaxz_spot = round(bbMax_z_spot[j]);

            ///test si spot appartient à cell
            if ( (bMinx_spot>=bMinx_cell) && (bMaxx_spot<=bMaxx_cell) && (bMiny_spot>=bMiny_cell) && (bMaxy_spot<=bMaxy_cell) && (bMinz_spot>=bMinz_cell) && (bMaxz_spot<=bMaxz_cell))
            {
                cpt++;
            }
        }

            int Rv;
            if(cpt>0 && cpt<=5)
            {
               Rv=2;
               genere_rayonfix_cell(cellSfr, cObj_cell, rBox_cell, Rv, N3D, cpt);
            }
            else if(cpt>5 && cpt<=10)
            {
               Rv=4;
               genere_rayonfix_cell(cellSfr, cObj_cell, rBox_cell, Rv, N3D, cpt);
            }
            else if(cpt>10 && cpt<=15)
            {
               Rv=6;
               genere_rayonfix_cell(cellSfr, cObj_cell, rBox_cell, Rv, N3D, cpt);
            }
            else if(cpt>15 && cpt<=20)
            {
               Rv=8;
               genere_rayonfix_cell(cellSfr, cObj_cell, rBox_cell, Rv, N3D, cpt);
            }
            else if(cpt>20 && cpt<=25)
            {
               Rv=10;
               genere_rayonfix_cell(cellSfr, cObj_cell, rBox_cell, Rv, N3D, cpt);
            }
            else if(cpt>25 && cpt<=30)
            {
               Rv=12;
               genere_rayonfix_cell(cellSfr, cObj_cell, rBox_cell, Rv, N3D, cpt);
            }
            else if(cpt>30)
            {
               Rv=14;
               genere_rayonfix_cell(cellSfr, cObj_cell, rBox_cell, Rv, N3D, cpt);
            }
    }

    ///calcul le centre de sphéroide
    int Cx_sphere = 0;
    int Cy_sphere = 0;
    int Cz_sphere = 0;
    Var3D cSphere = {Cx_sphere, Cy_sphere, Cz_sphere};
    Cx_sphere = (Cx_max - Cx_min)/2 + Cx_min;
    Cy_sphere = (Cy_max - Cy_min)/2 + Cy_min;
    Cz_sphere = (Cz_max - Cz_min);
    printf("Cx_sphere=%d,Cy_sphere=%d,Cz_sphere=%d\n", Cx_sphere, Cy_sphere, Cz_sphere);
    int Rx_sphere = 0;
    int Ry_sphere = 0;
    int R_sphere = 0;
    Rx_sphere = (Cx_max - Cx_min)/2;
    Ry_sphere = (Cy_max - Cy_min)/2;
    R_sphere = maximum(Rx_sphere, Ry_sphere, 0);
    printf("Rx_sphere=%d,Ry_sphere=%d, R_sphere=%d\n", Rx_sphere, Ry_sphere, R_sphere);

    ///export la carte des cellules en fonction de nombres de spots présents dans la cellule
    ///enregistrement binaire
    //std::string output_sfr = "E:\\" + userName + "\\" + acquiDate + "\\codeblocks\\cellSfr_" + oss.str() + ".bin";
    //Write_bin_ushort(cellSfr, N3D.x*N3D.y*N3D.z, output_sfr, "wb");
    ///enregistrement tiff
    std::string output_sfr = "E:\\" + userName + "\\" + acquiDate + "\\codeblocks\\cellSfr_" + oss.str() + ".tif";
    WRITE_tiff3D(cellSfr, output_sfr, N3D);

    ///calcul la distance euclidienne de chaque cellule par rapport au centre de sphéroide
    int dist_euc = 0;
    for(int i = 1; i <= lineNumber_cell; i ++){
        int Cx_cell = round(centerObject_x_cell[i]);
        int Cy_cell = round(centerObject_y_cell[i]);
        int Cz_cell = round(centerObject_z_cell[i]);
        int Rx_cell = round(rayon_x_cell[i]);
        int Ry_cell = round(rayon_y_cell[i]);
        int Rz_cell = round(rayon_z_cell[i]);
        int bMinx_cell = round(bbMin_x_cell[i]);
        int bMiny_cell = round(bbMin_y_cell[i]);
        int bMinz_cell = round(bbMin_z_cell[i]);
        int bMaxx_cell = round(bbMax_x_cell[i]);
        int bMaxy_cell = round(bbMax_y_cell[i]);
        int bMaxz_cell = round(bbMax_z_cell[i]);
        int rbox_x_cell = (bMaxx_cell -bMinx_cell)/2;
        int rbox_y_cell = (bMaxy_cell -bMiny_cell)/2;
        int rbox_z_cell = (bMaxz_cell -bMinz_cell)/2;
        int Rbox_max_cell = maximum(rbox_x_cell, rbox_y_cell, rbox_z_cell);
        int mInt_cell = round(meanIntensity_cell[i]);
        Var3D cObj_cell = {Cx_cell, Cy_cell, Cz_cell};
        Var3D rObj_cell = {Rx_cell, Ry_cell, Rz_cell};
        Var3D rBox_cell = {rbox_x_cell, rbox_y_cell, rbox_z_cell};

        dist_euc = round(sqrt((Cx_cell-Cx_sphere)*(Cx_cell-Cx_sphere)+(Cy_cell-Cy_sphere)*(Cy_cell-Cy_sphere)+(Cz_cell-Cz_sphere)*(Cz_cell-Cz_sphere)));
        genere_SFR_cell(distSfr, cObj_cell, rBox_cell, Rbox_max_cell, mInt_cell, N3D, dist_euc);

    }
    ///enregistrement binaire
    //std::string output_dist_euc = "E:\\" + userName + "\\" + acquiDate + "\\codeblocks\\dist_euc_" + oss.str() + ".bin";
    //Write_bin_ushort(distSfr, N3D.x*N3D.y*N3D.z, output_dist_euc, "wb");

    ///enregistrement tiff
    //std::string output_dist_euc = "E:\\" + userName + "\\" + acquiDate + "\\codeblocks\\dist_euc_" + oss.str() + ".tif";
    //WRITE_tiff3D(distSfr, output_dist_euc, N3D);

    ///calcul le nombre de spot
    int const nombrespot(10);
    int spot[nombrespot];
    int cpt0 = 0, cpt1 = 0, cpt2 = 0, cpt3 = 0, cpt4 = 0, cpt5 = 0, cpt6 = 0, cpt7 = 0, cpt8 = 0, cpt9 = 0;

    for(int i = 1; i <= lineNumber_spot; i ++){
        int Cx_spot = round(centerObject_x_spot[i]);
        int Cy_spot = round(centerObject_y_spot[i]);
        int Cz_spot = round(centerObject_z_spot[i]);

        int r_euc = 0;
        r_euc = sqrt((Cx_spot-Cx_sphere)*(Cx_spot-Cx_sphere)+(Cy_spot-Cy_sphere)*(Cy_spot-Cy_sphere)+(Cz_spot-Cz_sphere)*(Cz_spot-Cz_sphere));

        if ((r_euc > 0) && (r_euc <= 0.1*R_sphere))
        {
           cpt0++;
        }
        if ((r_euc > 0.1*R_sphere) && (r_euc <= 0.2*R_sphere))
        {
           cpt1++;
        }
        if ((r_euc > 0.2*R_sphere) && (r_euc <= 0.3*R_sphere))
        {
           cpt2++;
        }
        if ((r_euc > 0.3*R_sphere) && (r_euc <= 0.4*R_sphere))
        {
           cpt3++;
        }
        if ((r_euc > 0.4*R_sphere) && (r_euc <= 0.5*R_sphere))
        {
           cpt4++;
        }
        if ((r_euc > 0.5*R_sphere) && (r_euc <= 0.6*R_sphere))
        {
           cpt5++;
        }
        if ((r_euc > 0.6*R_sphere) && (r_euc <= 0.7*R_sphere))
        {
           cpt6++;
        }
        if ((r_euc > 0.7*R_sphere) && (r_euc <= 0.8*R_sphere))
        {
           cpt7++;
        }
        if ((r_euc > 0.8*R_sphere) && (r_euc <= 0.9*R_sphere))
        {
           cpt8++;
        }
        if ((r_euc > 0.9*R_sphere) && (r_euc <= R_sphere))
        {
           cpt9++;
        }
    }
    spot[9] = cpt9;
    spot[8] = cpt8;
    spot[7] = cpt7;
    spot[6] = cpt6;
    spot[5] = cpt5;
    spot[4] = cpt4;
    spot[3] = cpt3;
    spot[2] = cpt2;
    spot[1] = cpt1;
    spot[0] = cpt0;
    int somme = spot[9]+spot[8]+spot[7]+spot[6]+spot[5]+spot[4]+spot[3]+spot[2]+spot[1]+spot[0];
    std::cout<<"spot[9]="<<spot[9]<<std::endl;
    std::cout<<"spot[8]="<<spot[8]<<std::endl;
    std::cout<<"spot[7]="<<spot[7]<<std::endl;
    std::cout<<"spot[6]="<<spot[6]<<std::endl;
    std::cout<<"spot[5]="<<spot[5]<<std::endl;
    std::cout<<"spot[4]="<<spot[4]<<std::endl;
    std::cout<<"spot[3]="<<spot[3]<<std::endl;
    std::cout<<"spot[2]="<<spot[2]<<std::endl;
    std::cout<<"spot[1]="<<spot[1]<<std::endl;
    std::cout<<"spot[0]="<<spot[0]<<std::endl;
    std::cout<<"somme="<<somme<<std::endl;

    hist_spot << " T = " << t << std::endl;
    for(int i=0; i<10; i++)
    {
        hist_spot << spot[i] << std::endl;
    }

///fin de comment

    delete[] cellSfr;
    delete[] distSfr;
    }

    hist_spot.close();
    nombre_cell.close();
    nombre_spot.close();

}
