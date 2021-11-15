#include <cstdlib>
#include <iostream>

#define _USE_MATH_DEFINES
#include <math.h>
#include <cmath>

#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include <time.h>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image/stb_image_write.h"

////ASSIGNMENT 5- Filters- DIGITAL IMAGE - ANUSHA SHANKER///////////

using namespace std;

int mod(int a, int b) { return (a % b + b) % b; }

float vEkernel[3][3] = { {-1,0,1},{-1,0,1},{-1,0,1} };
float hEkernel[3][3] = { {1,1,1},{0,0,0},{-1,-1,-1} };
float MBkernel[5][5];

float edge[3][3] = { {1,0,-1},{0,0,0},{-1,0,1} };

int gridfactor = 1;
int T = 50;

float edge1[3][3] = { {0,1,0},{1,-4,1},{0,1,0} };
int width, height, channels1, width2, height2, channels2;
unsigned char* img2 = stbi_load("city.jpg", &width, &height, &channels2, 0);

unsigned char* diff;
unsigned char* diffp;
unsigned char* refI;
unsigned char* canvas;
unsigned char* dx;
unsigned char* dy;

unsigned char* gradx;
unsigned char* grady;

int minSL = 4;
int maxSL = 16;

float fc = 1;

struct Stroke
{
	int x, y;

};

vector<int> Kx;
vector<int>Ky;
vector<vector<int>> Sx;
vector<vector<int>> Sy;
vector<int>counter;

//Stroke S;


void genKernel(float MBkernel[5][5], int h, int w,int radius)
{
	//float MBkernel[3][3];
	//float sigma = 30;
	float sum = 0;


	//random filter kernel generation
	for (int i = 0; i < h; i++)
	{
		for (int j = 0; j < w; j++)

		{
			float r = sqrt(i * i + j * j);
			//float rad = 30;
			MBkernel[i][j] = exp(-(r * r) / (2 * radius * radius)) / (2 * M_PI) * radius * radius;
			if (MBkernel[i][j] < 0) { MBkernel[i][j] = 0; }
			sum += MBkernel[i][j];
		}
	}

	//normalizing
	for (int i = 0; i < h; i++)
	{
		for (int j = 0; j < w; j++)
		{

			MBkernel[i][j] = MBkernel[i][j] / sum;
		}
	}

	//cout << sum << endl;
	//checking
	float sum1 = 0;
	for (int i = 0; i < h; i++)
	{
		for (int j = 0; j < w; j++)
		{
			sum1 += MBkernel[i][j];
		}
	}
	//cout << sum1 << endl;
}

const int h = 5, w = 5, m = 3, n = 3;
int llim = -(w - 1) / 2;
int ulim = (w - 1) / 2;
int llim1 = -(m - 1) / 2;
int ulim1 = (n - 1) / 2;

int frowc = (w - 1) / 2;
int fcolc = (h - 1) / 2;
int frowc1 = (m - 1) / 2;
int fcolc1 = (n - 1) / 2;

void Guass(unsigned char* ref, unsigned char* img2, float MBkernel[5][5], int height, int width)
{
	
	for (int i = 0; i < height; i++)
	{


		for (int j = 0; j < width; j++)
		{
			int  Rval = 0, Gval = 0, Bval = 0;
			
			
			int rpos = (i * width + j) * 3;
			int gpos = rpos + 1;
			int bpos = rpos + 2;

			
				for (int a = llim; a <= ulim; a++)
				{
					for (int b = llim; b <= ulim; b++)
					{

						int x = j + b;
						int y = i + a;


						if (x<0 || x>width)
						{
							x = width - mod(a, width) - 1;

						}


						if (y<0 || y>height)
						{
							y = height - mod(b, height) - 1;
						}

						int rpos1 = (y * width + x) * 3;
						int gpos1 = rpos1 + 1;
						int bpos1 = rpos1 + 2;


						//cout << rpos << endl;
						Rval += img2[rpos1] * MBkernel[frowc + a][frowc + b];
						Gval += img2[gpos1] * MBkernel[frowc + a][frowc + b];
						Bval += img2[bpos1] * MBkernel[frowc + a][frowc + b];


					}

				




				
				}

			ref[rpos] = Rval;
			ref[gpos] = Gval;
			ref[bpos] = Bval;

		}
	}
}


void gradientx(unsigned char* gradx, unsigned char* ref, unsigned char* dx, float hEkernel[3][3], int height, int width)
{


	
	for (int i = 0; i < height; i++)
	{


		for (int j = 0; j < width; j++)
		{
			int  Rval = 0, Gval = 0, Bval = 0;

			int rpos = (i * width + j) * 3;
			int gpos = rpos + 1;
			int bpos = rpos + 2;
			int ind = (i * width + j);

			for (int a = llim1; a <= ulim1; a++)
			{
				for (int b = llim1; b <= ulim1; b++)
				{

					int x = j + b;
					int y = i + a;


					if (x<0 || x>width)
					{
						x = width - mod(a, width) - 1;

					}


					if (y<0 || y>height)
					{
						y = height - mod(b, height) - 1;
					}

					int rpos1 = (y * width + x) * 3;
					int gpos1 = rpos1 + 1;
					int bpos1 = rpos1 + 2;


					//cout << rpos << endl;
					Rval += ref[rpos1] * hEkernel[frowc1 + a][frowc1 + b];
					Gval += ref[gpos1] * hEkernel[frowc1 + a][frowc1 + b];
					Bval += ref[bpos1] * hEkernel[frowc1 + a][frowc1 + b];


				}

			}


			gradx[rpos] = (Rval + 255) / 2;
			gradx[gpos] = (Gval + 255) / 2;
			gradx[bpos] = (Bval + 255) / 2;

			dx[ind] = sqrt(pow(gradx[rpos], 2) + pow(gradx[gpos], 2) + pow(gradx[bpos], 2));
		}
	}



}

void gradienty(unsigned char* grady, unsigned char* ref, unsigned char* dy, float vEkernel[3][3], int height, int width)
{



	for (int i = 0; i < height; i++)
	{


		for (int j = 0; j < width; j++)
		{
			int  Rval = 0, Gval = 0, Bval = 0;

			int rpos = (i * width + j) * 3;
			int gpos = rpos + 1;
			int bpos = rpos + 2;

			int ind = (i * width + j);

			for (int a = llim1; a <= ulim1; a++)
			{
				for (int b = llim1; b <= ulim1; b++)
				{

					int x = j + b;
					int y = i + a;


					if (x<0 || x>width)
					{
						x = width - mod(a, width) - 1;

					}


					if (y<0 || y>height)
					{
						y = height - mod(b, height) - 1;
					}

					int rpos1 = (y * width + x) * 3;
					int gpos1 = rpos1 + 1;
					int bpos1 = rpos1 + 2;


					//cout << rpos << endl;
					Rval += ref[rpos1] * vEkernel[frowc1 + a][frowc1 + b];
					Gval += ref[gpos1] * vEkernel[frowc1 + a][frowc1 + b];
					Bval += ref[bpos1] * vEkernel[frowc1 + a][frowc1 + b];


				}

			}


			grady[rpos] = (Rval + 255) / 2;
			grady[gpos] = (Gval + 255) / 2;
			grady[bpos] = (Bval + 255) / 2;


			dy[ind] = sqrt(pow(gradx[rpos], 2) + pow(gradx[gpos], 2) + pow(gradx[bpos], 2));
		}
	}



}
void makeStroke(int radius, int x0, int y0, unsigned char* ref)
{

	int rpos = (x0 * width + y0) * 3;
	int gpos = rpos + 1;
	int bpos = rpos + 2;

	int strokR = ref[rpos];
	int strokG = ref[gpos];
	int strokB = ref[bpos];

	
	Kx.push_back(x0);
	Ky.push_back(y0);
	


	int x = x0;
	int y = y0;
	int lastdx = 0;
	int lastdy = 0;

	int ind = (x * width + y);

	int rpos1 = (x * width + y) * 3;
	int gpos1 = rpos1 + 1;
	int bpos1 = rpos1 + 2;

	diffp[ind] = sqrt(pow((ref[rpos1] - strokR), 2) + pow((ref[gpos1] - strokG), 2) + pow((ref[bpos1] - strokB), 2));



	for (int i = 1; i <= maxSL; i++)
	{
		if (i > minSL && abs(diff[ind]) < abs(diffp[ind]))
		{
			//cout << "check1" << endl;
			return;
		}

		if (dx[ind] == 0 && dy[ind] == 0)
		{
			//cout << "check2" << endl;
			return;
		}

		float gx = dx[ind];// (gradx[rpos1] + gradx[gpos1] + gradx[bpos1]) / 3;
		float gy = dy[ind];// (grady[rpos1] + grady[gpos1] + grady[bpos1]) / 3;
		float ddx = gx;
		float ddy = -gy;

		if (lastdx * ddx + lastdy * ddy < 0)
		{
			ddx = -ddx;
			ddy = -ddy;

		}

		ddx = fc * ddx + (1 - fc) * (lastdx);
		ddy = fc * ddy + (1 - fc) * (lastdy);

		float dxdy = sqrt(pow(ddx, 2) + pow(ddy, 2));

		ddx = ddx / dxdy;
		ddy = ddy / dxdy;

		x = x + radius * ddx;
		y = y + radius * ddy;

		lastdx = ddx;
		lastdy = ddy;

		if (x < 0 || y < 0 || x >= height || y >= width)
			break;
		else 
		{
			Kx.push_back(x);
			Ky.push_back(y);
		}


		

	}


	
	
	return;
	
	

}


void paintLayer(unsigned char* canvas, unsigned char* ref, int radius)
{
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			
			int rpos = (i * width + j) * 3;
			int gpos = rpos + 1;
			int bpos = rpos + 2;
			int ind = (i * width + j);

			diff[ind] = sqrt(pow((ref[rpos] - canvas[rpos]), 2) + pow((ref[gpos] - canvas[gpos]), 2) + pow((ref[bpos] - canvas[bpos]), 2));
		}
	}
	int grid = gridfactor * radius;
	int ct = 0;
	for (int x = 0; x < height; x += grid)
	{
		for (int y = 0; y < width; y += grid)
		{
			float areaErr = 0;
			int Maxx, Maxy, max = -333;
			
			for (int m = -grid / 2; m <= grid / 2; m++)
			{
				for (int n = -grid / 2; n <= grid / 2; n++)
				{
					int a = m + x;
					int b = n + y;

					int ind = (a * width + b);
					if (a >= 0 && b >=0 && a < height && b < width)
					{
						areaErr += diff[ind] / pow(grid, 2);

						if (diff[ind] > max)
						{
							Maxx = a;
							Maxy = b;
							int ind2 = (Maxx * width + Maxy);
							max = diff[ind2];
						}
					}


					


				}
			}

			if (areaErr > T)
			{

				(makeStroke(radius, Maxx, Maxy, ref));
				Sx.push_back(Kx);
				Sy.push_back(Ky);
				counter.push_back(ct++);

			}
			Kx.clear();
			Ky.clear();
			
			
			 
			//cout << S[x].y << endl;
			//cout << "testttt" << endl;
		}
		
	}
	
	cout << Sx.size() << endl;

	random_shuffle(counter.begin(), counter.end());


	for (int s = 0; s < counter.size(); s++)
	{
		for (int p = 0; p < Sx.at(counter.at(s)).size(); p++)
		{
			int xp = Sx.at(counter.at(s)).at(p);
			int yp = Sy.at(counter.at(s)).at(p); 

			int rp1 = (xp * width + yp) * 3;
			int gp1 = rp1 + 1;
			int bp1 = rp1 + 2;


			//cout << "tessssst" << endl;
			for (int j = 0; j < height; j++)
			{
				for (int i = 0; i < width; i++)
				{
					if ((pow((i - yp), 2) + pow((j - xp), 2) - pow(radius, 2)) <= 0)
					{
						int rp2 = (j * width + i) * 3;
						int gp2 = rp2 + 1;
						int bp2 = rp2 + 2;

						canvas[rp2] = ref[rp1];
						canvas[gp2] = ref[gp1];
						canvas[bp2] = ref[bp1];
					}

					
				}

			}
			//cout << "ugh" << endl;
			stbi_write_jpg("result9.jpg", width, height, 3, canvas, 100);
		}
	}
	
	Sx.clear();
	Sy.clear();
	
	//stbi_write_jpg("result9.jpg", width, height, 3, canvas, 100);
}








int main()
{
	
	

	

	//unsigned char* img1 = stbi_load("cat.jpg", &width, &height, &channels1, 0);
	
	//cout << "image loaded" << endl;
	diffp = new unsigned char[height * width];

	refI = new unsigned char[height * width * 3];
	canvas = new unsigned char[height * width * 3];
	dx = new unsigned char[height * width];
	dy = new unsigned char[height * width];
	gradx = new unsigned char[height * width * 3];
	grady = new unsigned char[height * width * 3];
	diff = new unsigned char[height * width];
	int radius = 16;
	int brushr[3] = { 8,4,2 };
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			int rpos = (i * width + j) * 3;
			int gpos = rpos + 1;
			int bpos = rpos + 2;
			

			canvas[rpos] = 255;
			canvas[gpos] = 255;
			canvas[bpos] = 255;

		}
	}



	

	for (int i = 0; i < 3; i++) 
	{

		genKernel(MBkernel, h, w, brushr[i]);
		

		Guass(refI, img2, MBkernel, height, width);
		gradientx(gradx, refI, dx, hEkernel, height, width);
		
		gradienty(grady, refI, dy, vEkernel, height, width);
		


		paintLayer(canvas, refI, brushr[i]);

		cout << "test4" << endl;
	}
	



















	
	//stbi_write_jpg("result2.jpg", width, height, 3, gradx, 100);
	//stbi_write_jpg("result3.jpg", width, height, 3, grady, 100);


}