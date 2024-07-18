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

typedef struct {
  int x,y;
}Var2D;

typedef struct {
  int x,y,z;
}Var3D;

class CSVRow
{
    public:
        std::string const& operator[](std::size_t index) const
        {
            return m_data[index];
        }
        std::size_t size() const
        {
            return m_data.size();
        }
        void readNextRow(std::istream& str)
        {
            std::string         line;
            std::getline(str, line);

            std::stringstream   lineStream(line);
            std::string         cell;

            m_data.clear();
            while(std::getline(lineStream, cell, ','))
            {
                m_data.push_back(cell);
            }
            // This checks for a trailing comma with no data after it.
            if (!lineStream && cell.empty())
            {
                // If there was a trailing comma then add an empty element.
                m_data.push_back("");
            }
        }
    private:
        std::vector<std::string>    m_data;
};

void Read_bin(std::string chemin, unsigned char resultat[], short precision, int NbPix);
void Write_bin_uchar(unsigned char *var_sav, int NbPix2D, std::string chemin, char options[]);
void Write_bin_ushort(unsigned short *var_sav, int NbPix2D, std::string chemin, char options[]);
void WRITE_tiff(unsigned short *var_sav, std::string chemin, int dim);
void WRITE_tiff3D(unsigned short *var_sav, std::string chemin, Var3D dim);
int maximum(int a, int b, int c);
void genere_rayonfix_cell(unsigned short *objet, Var3D cObj, Var3D rBox, int Rfixe, Var3D dim, int i);
void genere_SFR_cell(unsigned short *objet, Var3D cObj, Var3D rBox, int Rbox_max, int mInt, Var3D dim, int i);
void genere_BBox_cell(unsigned short *objet, Var3D bMin, Var3D bMax, Var3D cObj, int mInt, int sObj, Var3D dim, int i);
void genere_SFR_spot(unsigned short *objet, Var3D cObj, Var3D rBox, int Rbox_max, int mInt, Var3D dim, int i);
void genere_BBox_spot(unsigned short *objet, Var3D bMin, Var3D bMax, Var3D cObj, int mInt, int sObj, Var3D dim, int i);
void lire_csv(std::string userName, std::string acquiDate, int frame, Var3D N3D);

