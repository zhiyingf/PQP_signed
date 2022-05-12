/////////////////////////////////////////////////////////////////////////////////////////////
//	FileName:	MarchingCubes.cpp
//	Author	:	Michael Y. Polyakov
//	email	:	myp@andrew.cmu.edu  or  mikepolyakov@hotmail.com
//	website	:	www.angelfire.com/linux/myp
//	date	:	July 2002
//	
//	Description:	Marching Cubes Algorithm
/////////////////////////////////////////////////////////////////////////////////////////////

#pragma once
#include "MarchingCubes.h"
#include<fstream>
//#include <igl/readOBJ.h>
//#include <igl/readPLY.h>
//#include <imgui/imgui.h>

using namespace std;

mpVector LinearInterp(mp4Vector p1, mp4Vector p2, float value)
{
	mpVector p;
	if(p1.val != p2.val)
		p = (mpVector)p1 + ((mpVector)p2 - (mpVector)p1)/(p2.val - p1.val)*(value - p1.val);
	else 
		p = (mpVector)p1;
	return p;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//	MARCHING CUBES	//

//  VERSION  1A).  //
TRIANGLE* MarchingCubes(int ncellsX, int ncellsY, int ncellsZ, float minValue, mp4Vector * points,  
										INTERSECTION intersection, int &numTriangles)
{
	int xx = 3 * ncellsX * ncellsY * ncellsZ;
	TRIANGLE * triangles = new TRIANGLE[xx];//this should be enough space, if not change 4 to 5
	/*F = Eigen::MatrixXf::Zero(xx, 3);*/
	ofstream out("spereOneMC.obj");
	//out << "#mcTriangle.obj" << endl;

	numTriangles = int(0);
	int countpoint = 1;

	int YtimeZ = (ncellsY+1)*(ncellsZ+1);
	//go through all the points
	for(int i=0; i < ncellsX; i++)			//x axis
		for(int j=0; j < ncellsY; j++)		//y axis
			for(int k=0; k < ncellsZ; k++)	//z axis
			{
				//initialize vertices
				mp4Vector verts[8];
				int ind = i*YtimeZ + j*(ncellsZ+1) + k;
				//int ind = k * (ncellsY + 1) * (ncellsX + 1) + j * (ncellsX + 1) + i;
   /*(step 3)*/ verts[0] = points[ind];//(0,0,0)
				verts[1] = points[ind + YtimeZ];//(1,0,0)
				verts[2] = points[ind + YtimeZ + 1];//(1,0,1)
				verts[3] = points[ind + 1];//(0,0,1)
				verts[4] = points[ind + (ncellsZ+1)];//(0,1,0)
				verts[5] = points[ind + YtimeZ + (ncellsZ+1)];//(1,1,0)
				verts[6] = points[ind + YtimeZ + (ncellsZ+1) + 1];//(1,1,1)
				verts[7] = points[ind + (ncellsZ+1) + 1];//(0,1,1)

				////int ind = i * YtimeZ + j * (ncellsZ + 1) + k;
				//int ind = k * (ncellsY + 1) * (ncellsX + 1) + j * (ncellsX + 1) + i;
				///*(step 3)*/ verts[0] = points[ind];
				///*verts[1] = points[ind + YtimeZ];
				//verts[2] = points[ind + YtimeZ + (ncellsZ + 1)];
				//verts[3] = points[ind + (ncellsZ + 1)];
				//verts[4] = points[ind + 1] ;
				//verts[5] = points[ind + YtimeZ + 1] ;
				//verts[6] = points[ind + YtimeZ + (ncellsZ + 1) + 1];
				//verts[7] = points[ind + (ncellsZ + 1) + 1];*/

				//	verts[1] = points[ind + YtimeZ];//(1,0,0)
				//verts[2] = points[ind + YtimeZ + 1];//(1,0,1)
				//verts[3] = points[ind + 1];//(0,0,1)
				//verts[4] = points[ind + (ncellsZ+1)];//(0,1,0)
				//verts[5] = points[ind + YtimeZ + (ncellsZ+1)];//(1,1,0)
				//verts[6] = points[ind + YtimeZ + (ncellsZ+1) + 1];//(1,1,1)
				//verts[7] = points[ind + (ncellsZ+1) + 1];//(0,1,1)


					
				
				//get the index
				int cubeIndex = int(0);
				for(int n=0; n < 8; n++)
   /*(step 4)*/		if(verts[n].val < minValue) cubeIndex |= (1 << n);

				//check if its completely inside or outside
   /*(step 5)*/ if(!edgeTable[cubeIndex]) continue;
			
				//get intersection vertices on edges and save into the array
   				mpVector intVerts[12];
				int indexpoint[12];
				
				/*(step 6)*/ if (edgeTable[cubeIndex] & 1) {
					intVerts[0] = intersection(verts[0], verts[1], minValue);
					indexpoint[0] = countpoint++;
					out << "v " << intVerts[0].x << " " << intVerts[0].y << " " << intVerts[0].z << endl;
                }
				if (edgeTable[cubeIndex] & 2) {
					intVerts[1] = intersection(verts[1], verts[2], minValue);
					indexpoint[1] = countpoint++;
					out << "v " << intVerts[1].x << " " << intVerts[1].y << " " << intVerts[1].z << endl;
				}
				if (edgeTable[cubeIndex] & 4) {
					intVerts[2] = intersection(verts[2], verts[3], minValue);
					indexpoint[2] = countpoint++;
					out << "v " << intVerts[2].x << " " << intVerts[2].y << " " << intVerts[2].z << endl;
				}
				if (edgeTable[cubeIndex] & 8) {
					intVerts[3] = intersection(verts[3], verts[0], minValue);
					indexpoint[3] = countpoint++;
					out << "v " << intVerts[3].x << " " << intVerts[3].y << " " << intVerts[3].z << endl;
				}
				if (edgeTable[cubeIndex] & 16) {
					intVerts[4] = intersection(verts[4], verts[5], minValue);
					indexpoint[4] = countpoint++;
					out << "v " << intVerts[4].x << " " << intVerts[4].y << " " << intVerts[4].z << endl;
				}
				if (edgeTable[cubeIndex] & 32) {
					intVerts[5] = intersection(verts[5], verts[6], minValue);
					indexpoint[5] = countpoint++;
					out << "v " << intVerts[5].x << " " << intVerts[5].y << " " << intVerts[5].z << endl;
				}
				if (edgeTable[cubeIndex] & 64) {
					intVerts[6] = intersection(verts[6], verts[7], minValue);
					indexpoint[6] = countpoint++;
					out << "v " << intVerts[6].x << " " << intVerts[6].y << " " << intVerts[6].z << endl;
				}
				if (edgeTable[cubeIndex] & 128) {
					intVerts[7] = intersection(verts[7], verts[4], minValue);
					indexpoint[7] = countpoint++;
					out << "v " << intVerts[7].x << " " << intVerts[7].y << " " << intVerts[7].z << endl;
				}
				if (edgeTable[cubeIndex] & 256) {
					intVerts[8] = intersection(verts[0], verts[4], minValue);
					indexpoint[8] = countpoint++;
					out << "v " << intVerts[8].x << " " << intVerts[8].y << " " << intVerts[8].z << endl;
				}
				if (edgeTable[cubeIndex] & 512) {
					intVerts[9] = intersection(verts[1], verts[5], minValue);
					indexpoint[9] = countpoint++;
					out << "v " << intVerts[9].x << " " << intVerts[9].y << " " << intVerts[9].z << endl;
				}
				if (edgeTable[cubeIndex] & 1024) {
					intVerts[10] = intersection(verts[2], verts[6], minValue);
					indexpoint[10] = countpoint++;
					out << "v " << intVerts[10].x << " " << intVerts[10].y << " " << intVerts[10].z << endl;
				}
				if (edgeTable[cubeIndex] & 2048) {
					intVerts[11] = intersection(verts[3], verts[7], minValue);
					indexpoint[11] = countpoint++;
					out << "v " << intVerts[11].x << " " << intVerts[11].y << " " << intVerts[11].z << endl;
				}

				//now build the triangles using triTable
				for (int n=0; triTable[cubeIndex][n] != -1; n+=3) {
					    int p1 = triTable[cubeIndex][n + 2];
						int p2 = triTable[cubeIndex][n + 1];
						int p3 = triTable[cubeIndex][n];

   /*(step 7)*/ 		triangles[numTriangles].p[0] = intVerts[p1];
						triangles[numTriangles].p[1] = intVerts[p2];
						triangles[numTriangles].p[2] = intVerts[p3];

   /*(step 8)*/ 		triangles[numTriangles].norm = ((triangles[numTriangles].p[1] - 
						triangles[numTriangles].p[0]).Cross(triangles[numTriangles].p[2] - 
						triangles[numTriangles].p[0])).Normalize();
						out << "f " << indexpoint[p1] << " " << indexpoint[p2] << " " << indexpoint[p3] << endl;
						numTriangles++;
				}
			
			}	//END OF FOR LOOP
		        
		out.close();
        //free all the wasted space
		TRIANGLE * retTriangles = new TRIANGLE[numTriangles];
		for(int i=0; i < numTriangles; i++)
			retTriangles[i] = triangles[i];
		delete [] triangles;
	
	return retTriangles;
}


//	VERSION  1B).  //
TRIANGLE* MarchingCubesLinear(int ncellsX, int ncellsY, int ncellsZ, float minValue, 
									mp4Vector * points, int &numTriangles)
{
	return MarchingCubes(ncellsX, ncellsY, ncellsZ, minValue, points, LinearInterp, numTriangles);
}


//	VERSION  2A).  //
TRIANGLE* MarchingCubes(float mcMinX, float mcMaxX, float mcMinY, float mcMaxY, float mcMinZ, float mcMaxZ, 
							int ncellsX, int ncellsY, int ncellsZ, float minValue, 
							FORMULA formula, INTERSECTION intersection, int &numTriangles)
{
	//space is already defined and subdivided, staring with step 3
	//first initialize the points
	mp4Vector * mcDataPoints = new mp4Vector[(ncellsX+1)*(ncellsY+1)*(ncellsZ+1)];
	mpVector stepSize((mcMaxX-mcMinX)/ncellsX, (mcMaxY-mcMinY)/ncellsY, (mcMaxZ-mcMinZ)/ncellsZ);
	
	int YtimesZ = (ncellsY+1)*(ncellsZ+1);	//for extra speed
	for(int i=0; i < ncellsX+1; i++) {
		int ni = i*YtimesZ;						//for speed
		float vertX = mcMinX + i*stepSize.x;
		for(int j=0; j < ncellsY+1; j++) {
			int nj = j*(ncellsZ+1);				//for speed
			float vertY = mcMinY + j*stepSize.y;
			for(int k=0; k < ncellsZ+1; k++) {
				mp4Vector vert(vertX, vertY, mcMinZ + k*stepSize.z, 0);
				vert.val = formula((mpVector)vert);
   /*(step 3)*/ mcDataPoints[ni + nj + k] = vert;
			}
		}
	}
	//then run Marching Cubes (version 1A) on the data
	return MarchingCubes(ncellsX, ncellsY, ncellsZ, minValue, mcDataPoints, intersection, numTriangles);
}

//	VERSION  2B).  //
TRIANGLE* MarchingCubesLinear(float mcMinX, float mcMaxX, float mcMinY, float mcMaxY, float mcMinZ, float mcMaxZ, 
								int ncellsX, int ncellsY, int ncellsZ, float minValue, 
								FORMULA formula, int &numTriangles)
{
	return MarchingCubes(mcMinX, mcMaxX, mcMinY, mcMaxY, mcMinZ, mcMaxZ, ncellsX, ncellsY, ncellsZ, minValue,
		formula, LinearInterp, numTriangles);
}