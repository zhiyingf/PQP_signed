#include "BaseModel.h"
#include "Distance_OBB.h"
#include "MarchingCubes.h"
#include <fstream>

#include <io.h> 
#include <regex>
#include <filesystem>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/Point_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel k;
typedef k::Point_3 Point3;
typedef CGAL::Polyhedron_3<k> poly;
typedef CGAL::Polyhedral_mesh_domain_3<poly, k> mesh_domain;


const float minValue = 0;
TRIANGLE* Triangles;
int numOfTriangles;
mp4Vector* mcPoints;

//set bounding box
#define MINX -0.5
#define MAXX 0.5
#define MINY -0.5
#define MAXY 0.5
#define MINZ -0.5
#define MAXZ 0.5

using namespace std;
//namespace fs = std::filesystem;

const int nX = 100;
const int nY = 100;
const int nZ = 100;



PQP_REAL computeDis(PQP::Distance_OBB* obb, Model3D::CPoint3D pos, mesh_domain* domain)
{
    //使用cgal判断point在model内外
    //自行判断point在model内外：最近面法线与pp'的点乘
    Point3 p(pos.x, pos.y, pos.z);
    bool b = domain->is_in_domain_object()(p) > 0;

    PQP::Distance_OBB::QueryResult intersects;
    intersects = (*obb).Query(pos);
    PQP_REAL dis = intersects.distance;
    if (!b) dis = -dis;
    return dis;
}

void InitDataMc(PQP::Distance_OBB* obb, mesh_domain* domain, string fileName) {
    int size = (nX + 1) * (nY + 1) * (nZ + 1);
    mcPoints = new mp4Vector[size];


    float* sdf = new float[size];


    mpVector stepSize((MAXX - MINX) / nX, (MAXY - MINY) / nY, (MAXZ - MINZ) / nZ);
    for (int i = 0; i <= nX; i++)
        for (int j = 0; j <= nY; j++)
            for (int k = 0; k <= nZ; k++) {
                mp4Vector vert(MINX + i * stepSize.x, MINY + j * stepSize.y, MINZ + k * stepSize.z, 0);

                Model3D::CPoint3D pos(MINX + i * stepSize.x, MINY + j * stepSize.y, MINZ + k * stepSize.z);

                vert.val = computeDis(obb, pos, domain);

                mcPoints[i * (nY + 1) * (nZ + 1) + j * (nZ + 1) + k] = vert;

                //sdf[i * (nY + 1) * (nZ + 1) + j * (nZ + 1) + k] = vert.val;

                //(x,z)---(z,x)---(nX-1-x,z)互换x,z且（x,z）逆时针旋转90°
                int idx = (nX - i) + j * (nZ + 1) + k * (nY + 1) * (nZ + 1);//(x,z)---(z,x)---(nX-1-x,z)

                sdf[idx] = -vert.val;
            }


    ///float array to binary file
    string str1 = string("./sdf/") + fileName + "-100.txt";
    std::ofstream  ofs(str1, std::ios::binary | std::ios::out);
    ofs.write((const char*)sdf, sizeof(float) * size);
    ofs.close();

}

//get *.obj file
void getFiles(std::string path, std::vector<std::string>& files, std::vector<std::string>& filenames)
{
    intptr_t   hFile = 0;//intptr_t和uintptr_t的类型:typedef long int； typedef unsigned long int  
    struct _finddata_t fileinfo;
    std::string p;
    if ((hFile = _findfirst(p.assign(path).append("/*").c_str(), &fileinfo)) != -1)//assign方法：先将原字符串清空，然后赋予新的值作替换。  
    {
        do
        {
            if (strcmp(fileinfo.name, ".") != 0 && strcmp(fileinfo.name, "..") != 0)
            {
                files.push_back(p.assign(path).append("/").append(fileinfo.name));
                filenames.push_back(fileinfo.name);
            }
        } while (_findnext(hFile, &fileinfo) == 0);
        _findclose(hFile);
    }
}

//compute SDF - save SDF binary file - verify SDF(ie.extract surface using Marching cubes)
void solved(const vector<string>& fileVec, const vector<string>& fileNam) {
    int i = 0;
    for (auto& v : fileVec) {
        //模型初始化
        /*string fileName("F:\\xinCode\\PQP_signed\\data\\all\\");
        fileName += v + ".off";*/
        string fileName = v;
        Model3D::CBaseModel model(fileName);//*.obj / *.off / *.m
        model.LoadModel();
        PQP::Distance_OBB* obb = new PQP::Distance_OBB(model);

        poly pol;
        mesh_domain* domain;
        ifstream mod(fileName);
        mod >> pol;
        mod.close();
        domain = new mesh_domain(pol);

        string name_model = fileNam[i++];
        name_model = name_model.substr(0, name_model.size() - 4);
        InitDataMc(obb, domain, name_model);
        Triangles = MarchingCubes(nX, nY, nZ, minValue, mcPoints, LinearInterp, numOfTriangles);
        system("pause");
    }
}

int main()
{
    //vector<string> fileVec = { "0","4","26","84","85","88" };

    vector<string> fileVec;
    vector<string> fileNam;
    //需要读取的文件夹路径，使用单右斜杠“/”  
    string filePath = "F:/xinCode/PQP_signed/data";
    //string filePath = "./datas"; //相对路径也可以  
    getFiles(filePath, fileVec, fileNam);

    solved(fileVec, fileNam);
    return 0;
}