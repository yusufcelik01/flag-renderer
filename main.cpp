#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#define _USE_MATH_DEFINES
#include <math.h>
#include <GL/glew.h>
//#include <OpenGL/gl3.h>   // The GL Header File
#include <GLFW/glfw3.h> // The GLFW header
#include <glm/glm.hpp> // GL Math library header
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp> 
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>
#include <glm/gtx/normal.hpp>

#include "helpers.hpp"

#define BUFFER_OFFSET(i) ((char*)NULL + (i))

using namespace std;

GLuint gProgram[2];
int gWidth, gHeight;

GLint modelingMatrixLoc[2];
GLint viewingMatrixLoc[2];
GLint projectionMatrixLoc[2];
GLint eyePosLoc[2];

glm::mat4 projectionMatrix;
glm::mat4 viewingMatrix;
glm::mat4 modelingMatrix;
glm::vec3 eyePos(0, 0, 0);

int activeProgramIndex = 0;

struct Vertex
{
    Vertex() : x(0), y(0), z(0) {}
    Vertex(GLfloat inX, GLfloat inY, GLfloat inZ) : x(inX), y(inY), z(inZ) { }
    Vertex(glm::vec3 vec) : x(vec.x), y(vec.y), z(vec.z) { }
    glm::vec3 toGlmVec3(){
        glm::vec3 v;
        v.x = x;
        v.y = y;
        v.z = z;
        
        return v;
    }

    Vertex operator+(Vertex v)
    {
        Vertex temp;
        temp.x = v.x + this->x;
        temp.y = v.y + this->y;
        temp.z = v.z + this->z;

        return temp;
    }

    Vertex operator-(Vertex v)
    {
        Vertex temp;
        temp.x = this->x - v.x;
        temp.y = this->y - v.y;
        temp.z = this->z - v.z;

        return temp;
    }

    Vertex operator*(GLfloat s)
    {
        Vertex temp;
        temp.x = s * this->x;
        temp.y = s * this->y;
        temp.z = s * this->z;

        return temp;
    }

    Vertex mirrorOver(Vertex mirror)
    {
        return mirror + (mirror - *this);
    }

    GLfloat x, y, z;
};




struct Texture
{
    Texture() : u(0), v(0) { }
    Texture(GLfloat inU, GLfloat inV) : u(inU), v(inV) { }
    GLfloat u, v;
    
    Texture operator+(Texture t)
    {
        Texture temp;
        temp.u = this->u + t.u;
        temp.v = this->v + t.v;

        return temp;
    }

    Texture operator/(GLfloat s)
    {
        Texture temp;
        temp.u = this->u/s; 
        temp.v = this->v/s;

        return temp;
    }

    Texture operator*(GLfloat s)
    {
        Texture temp;
        temp.u = this->u*s; 
        temp.v = this->v*s;

        return temp;
    }
};

struct Normal
{
    GLfloat x, y, z;

    Normal() : x(0), y(0), z(0) {}
    Normal(GLfloat inX, GLfloat inY, GLfloat inZ) : x(inX), y(inY), z(inZ) { }
    //Normal(glm::vec3 vec) : x(vec.x), y(vec.y), z(vec.z) { }
    Normal(glm::vec3 vec)
    {
        glm::vec3 norm = glm::normalize(vec);
        x = norm.x;
        y = norm.y;
        z = norm.z;

    } 
    Normal operator=(glm::vec3 v)
    {
        glm::vec3 norm = glm::normalize(v);
        x = norm.x;
        y = norm.y;
        z = norm.z;

        return *this;
    }
};

struct Face
{
	Face(int v[], int t[], int n[]) {
		vIndex[0] = v[0];
		vIndex[1] = v[1];
		vIndex[2] = v[2];
		tIndex[0] = t[0];
		tIndex[1] = t[1];
		tIndex[2] = t[2];
		nIndex[0] = n[0];
		nIndex[1] = n[1];
		nIndex[2] = n[2];
	}
    GLuint vIndex[3], tIndex[3], nIndex[3];
};

struct BezierPatch
{
    Vertex controlPoints[4][4];
    Texture texTopLeft, texTopRight, texBottomRight, texBottomLeft;
};

vector<Vertex> gVertices;
vector<Texture> gTextures;
vector<Normal> gNormals;
vector<Face> gFaces;
vector<BezierPatch> gPatches;

GLuint gVertexAttribBuffer, gIndexBuffer;
GLint gInVertexLoc, gInNormalLoc;
int gVertexDataSizeInBytes, gNormalDataSizeInBytes;
int gTextureDataSizeInBytes;

char* imagePath = NULL;
unsigned char*  rawImage = NULL;
unsigned int flagTexture;
float flagAspecRat = 1.6;

int sampleSize = 10;
int patchSize = 1;

GLfloat* linspace(GLfloat start, GLfloat end, int n) 
{
    assert(n > 0); 
    //GLfloat* arr = (GLfloat*) calloc(n, sizeof(GLfloat));
    GLfloat* arr = new GLfloat [n];

    GLfloat step = (end-start)/(n-1);
    for(int i = 0; i<n; i++, start+=step)
    {
       arr[i] = start; 
    }

    return arr;
}

int factorial(int n)
{
    if(n < 2) {
        return 1;
    }

    return n * factorial(n-1);
}

int nCombinationK(int n, int k)
{
    //TODO optimize

    //if(k*2 > n) {
    //    k = n - k;
    //}   
    //int combination = 1;
    //while(k){
    //    combination *= n;
    //    combination /= k;
    //    n--;
    //    k--;
    //}
    
    //return combination;
    return factorial(n) / factorial(k) / factorial(n-k);
}

GLfloat bern(int i, int n, GLfloat t)
{
    //TODO optimize

    if(i > n)
    {
        return -1;
    }

    return nCombinationK(n, i) * pow(t, i) * pow(1-t, n-i);
}

bool sampleBezier()
{
    Vertex controlPoints[4][4];
    
    //controlPoints[0][0] = Vertex(0, 0,0);
    //controlPoints[0][1] = Vertex(0,-1,0);
    //controlPoints[0][2] = Vertex(0,-2,0);
    //controlPoints[0][3] = Vertex(0,-3,0);

    //controlPoints[1][0] = Vertex(3,0,1);
    //controlPoints[1][1] = Vertex(3,-1,2);
    //controlPoints[1][2] = Vertex(3,-2,2);
    //controlPoints[1][3] = Vertex(3,-3,1);

    //controlPoints[2][0] = Vertex(6,0,-1);
    //controlPoints[2][1] = Vertex(6,-1,-2);
    //controlPoints[2][2] = Vertex(6,-2,-2);
    //controlPoints[2][3] = Vertex(6,-3,-1);
    //
    //controlPoints[3][0] = Vertex(8,0,0);
    //controlPoints[3][1] = Vertex(8,-1,0);
    //controlPoints[3][2] = Vertex(8,-2,0);
    //controlPoints[3][3] = Vertex(8,-3,0);


    controlPoints[0][0] = Vertex(0, 0,0);
    controlPoints[0][1] = Vertex(0,-1,0);
    controlPoints[0][2] = Vertex(0,-2,0);
    controlPoints[0][3] = Vertex(0,-3,0);

    controlPoints[1][0] = Vertex(3,0,3);
    controlPoints[1][1] = Vertex(3,-1,5);
    controlPoints[1][2] = Vertex(3,-2,6);
    controlPoints[1][3] = Vertex(3,-3,4);

    controlPoints[2][0] = Vertex(6,0,-4);
    controlPoints[2][1] = Vertex(6,-1,-6);
    controlPoints[2][2] = Vertex(6,-2,-5);
    controlPoints[2][3] = Vertex(6,-3,-3);

    controlPoints[3][0] = Vertex(8,0,0);
    controlPoints[3][1] = Vertex(8,-1,0);
    controlPoints[3][2] = Vertex(8,-2,0);
    controlPoints[3][3] = Vertex(8,-3,0);


    Vertex **Q; 

    Q = new Vertex*[sampleSize];
    GLfloat* samplePoints;

    samplePoints = linspace(0, 1, sampleSize);
    for(int row = 0; row < sampleSize; row++)
    {
        Q[row] = new Vertex[sampleSize];
        //Qy[row] = new Vertex[sampleSize];
        //Qz[row] = new Vertex[sampleSize];

        memset(&Q[row][0], 0, sizeof(struct Vertex)*sampleSize);
        //TODO remove logging debug
        //for(int col = 0; col< sampleSize; col++)
        //{
        //std::cout << "Q " << row << " " << col<<" init: " <<" "<< Q[row][col].x << " "<< Q[row][col].y <<" "<< Q[row][col].z << std::endl;
        //}
    }

    for(int s = 0; s < sampleSize; s++){
        for(int t = 0; t < sampleSize; t++){
            
            for(int i=0; i<=3; i++){
                for(int j = 0; j<=3; j++){
                    Q[s][t].x += bern(i, 3, samplePoints[s]) * bern(j, 3, samplePoints[t]) * controlPoints[i][j].x;
                    Q[s][t].y += bern(i, 3, samplePoints[s]) * bern(j, 3, samplePoints[t]) * controlPoints[i][j].y;
                    Q[s][t].z += bern(i, 3, samplePoints[s]) * bern(j, 3, samplePoints[t]) * controlPoints[i][j].z;
                    //TODO optimize bern function is called 3 times
                    //bern function can be computed more effiecntly
                }
            }

        }
    }

    
    for(int s = 0; s < sampleSize; s++){
        for(int t = 0; t < sampleSize; t++){
            gVertices.push_back(Q[s][t]);
        }
    }

    vector<glm::vec3> faceNormals; 
   
    for(int s = 0; s <= sampleSize-2; s++)
    {
        for(int t = 0; t <= sampleSize-2; t++)
        {
            glm::vec3 v1, v2, v3, v4; 
            //glm::vec3 norm1, norm2;

            int vIndex[3], tIndex[3];

            vIndex[0] = s*sampleSize + t;
            vIndex[1] = s*sampleSize + t + 1;
            vIndex[2] = (s+1)*sampleSize + t;

            //v1 = glm::vec3(Q[s][t].x, Q[s][t].y, Q[s][t].z);
            //v2 = glm::vec3(Q[s][t+1].x, Q[s][t+1].y, Q[s][t+1].z);
            //v3 = glm::vec3(Q[s+1][t].x, Q[s+1][t].y, Q[s+1][t].z);
            //v4 = glm::vec3(Q[s+1][t+1].x, Q[s+1][t+1].y, Q[s+1][t+1].z);

            v1 = Q[s][t].toGlmVec3();
            v2 = Q[s][t+1].toGlmVec3();
            v3 = Q[s+1][t].toGlmVec3();
            v4 = Q[s+1][t+1].toGlmVec3();

            //norm1 = cross(v2-v1, v3-v2);//TODO optimize
            //norm2 = cross(v3-v2, v4-v2);//v3-v2 is computed twice
            
            //faceNormals.push_back(norm1);
            //faceNormals.push_back(norm2);
            faceNormals.push_back(glm::triangleNormal(v1,v2,v3));
            faceNormals.push_back(glm::triangleNormal(v3,v2,v4));

            int ind = faceNormals.size()-2;
            glm::vec3 tempv = faceNormals[ind];
            //std::cout << "face normal " << ind << "" << tempv.x << " " << tempv.y <<" "<< tempv.z << std::endl;
            ind++;
            tempv = faceNormals[ind];
            //std::cout << "face normal " << ind << "" << tempv.x << " " << tempv.y <<" "<< tempv.z << std::endl;

            gFaces.push_back(Face(vIndex, tIndex, vIndex));
            vIndex[0] = (s+1)*sampleSize + t;
            vIndex[1] = s*sampleSize + t + 1;
            vIndex[2] = (s+1)*sampleSize + t + 1;

            gFaces.push_back(Face(vIndex, vIndex, vIndex));
        }
    }


    gNormals.resize((sampleSize)*(sampleSize));

    gNormals[0] = Normal(faceNormals[0]);//top-left
    //std::cout << "normal " << faceNormals[(sampleSize-1)*2 -1].x << " "  << faceNormals[(sampleSize-1)*2 -1].y << " " << faceNormals[(sampleSize-1)*2 -1].z << " "<< std::endl;
    //std::cout << "normal " << faceNormals[(sampleSize-1 )*2 -2].x << " " << faceNormals[(sampleSize-1 )*2 -2].y << " " << faceNormals[(sampleSize-1 )*2 -2].z << " "<< std::endl;
    gNormals[sampleSize-1] = Normal((faceNormals[(sampleSize-1)*2 -1] + faceNormals[(sampleSize-1)*2 -2])/2.f);//bottom-left corner vertex

    gNormals[(sampleSize)*(sampleSize-1)] = Normal(
            (faceNormals[2*(sampleSize-2)*(sampleSize-1)] 
             + faceNormals[2*(sampleSize-2)*(sampleSize-1) + 1])/2.f
            );//top right

    gNormals[sampleSize*sampleSize -1] = Normal(
            faceNormals[2*(sampleSize-1) * (sampleSize-1)-1]);//bottom right

    //compute vertex normals at the edges
    //for(int i = 1; i <= sampleSize-2; i++)
    //{
    //    int temp;
    //    gNormals[i] = Normal( (faceNormals[i*2 - 2] +
    //                faceNormals[i*2 -1] + 
    //                faceNormals[i*2]) / 3.f);//left edge of the patch

    //    temp = ((sampleSize-1)*(sampleSize-2)+ i)*2;  
    //    gNormals[sampleSize*(sampleSize-1) -1 + i] = Normal(
    //               (faceNormals[temp - 1] + 
    //                faceNormals[temp] + 
    //                faceNormals[temp + 1]) /3.f  );
    //    std::cout << "Vnorm index: " << sampleSize*(sampleSize-1) -1 + i 
    //        << " " << temp -1 << " " << temp << " " << temp+1
    //        <<std::endl;


    //    gNormals[(sampleSize*i)] = Normal(faceNormals[sampleSize * i]);//TODO   
    //    //TODO something seems fishy
    //    gNormals[(sampleSize*i) + sampleSize -1] = Normal(faceNormals[sampleSize * i + sampleSize -1]);
    //}
    for(int i = 1; i <= sampleSize-2; i++){
        int temp;
        //left edge of the patch
        gNormals[i] = Normal( (faceNormals[i*2 - 2] +
                    faceNormals[i*2 -1] + 
                    faceNormals[i*2]) / 3.f);

        //right edge of the patch
        temp = (sampleSize-1)*2*(sampleSize-2) + 2*i;
        gNormals[sampleSize*(sampleSize-1) + i] = Normal(
                    (faceNormals[temp-1] +
                     faceNormals[temp] +
                     faceNormals[temp+1]
                    )/3.f
                );
        //top edge vertex normals
        gNormals[sampleSize*i] = Normal(
                    (faceNormals[(sampleSize-1)*(i-1)] +
                     faceNormals[(sampleSize-1)*(i-1) + 1] +
                     faceNormals[(sampleSize-1) * i]
                    )/3.f
                );

        //bottom edge vertex normals
        gNormals[sampleSize*(i+1) - 1] = Normal(
                    (faceNormals[(sampleSize-1)*2*i -1] +
                     faceNormals[(sampleSize-1)*2*(i+1) -1] +
                     faceNormals[(sampleSize-1)*2*(i+1) -2]
                    )/3.f
                );

    }

    for(int s = 1; s <= sampleSize-2; s++)
    {
        for(int t = 1; t <= sampleSize-2; t++)
        {
            int temp;
            glm::vec3 norm;

            temp = ((sampleSize-1)*s + t) *2;
            norm = faceNormals[temp - 1] +
                        faceNormals[temp] +
                        faceNormals[temp +1];
            norm /= 3;
            gNormals[sampleSize*s + t] = Normal(norm);
        }

    }

    


    //for(int patchNo = 0; patchNo < 1; patchNo++) 
    //{
    //    //gVertices.push_back(Vertex());

    //}


    delete[] samplePoints;
    for(int row = 0; row < sampleSize; row++)
    {
        delete[] Q[row];
    }

    return true;
}



bool ReadDataFromFile(
    const string& fileName, ///< [in]  Name of the shader file
    string&       data)     ///< [out] The contents of the file
{
    fstream myfile;

    // Open the input 
    myfile.open(fileName.c_str(), std::ios::in);

    if (myfile.is_open())
    {
        string curLine;

        while (getline(myfile, curLine))
        {
            data += curLine;
            if (!myfile.eof())
            {
                data += "\n";
            }
        }

        myfile.close();
    }
    else
    {
        return false;
    }

    return true;
}

GLuint createVS(const char* shaderName)
{
    string shaderSource;

    string filename(shaderName);
    if (!ReadDataFromFile(filename, shaderSource))
    {
        cout << "Cannot find file name: " + filename << endl;
        exit(-1);
    }

    GLint length = shaderSource.length();
    const GLchar* shader = (const GLchar*) shaderSource.c_str();

    GLuint vs = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vs, 1, &shader, &length);
    glCompileShader(vs);

    char output[1024] = {0};
    glGetShaderInfoLog(vs, 1024, &length, output);
    printf("VS compile log: %s\n", output);

	return vs;
}

GLuint createFS(const char* shaderName)
{
    string shaderSource;

    string filename(shaderName);
    if (!ReadDataFromFile(filename, shaderSource))
    {
        cout << "Cannot find file name: " + filename << endl;
        exit(-1);
    }

    GLint length = shaderSource.length();
    const GLchar* shader = (const GLchar*) shaderSource.c_str();

    GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fs, 1, &shader, &length);
    glCompileShader(fs);

    char output[1024] = {0};
    glGetShaderInfoLog(fs, 1024, &length, output);
    printf("FS compile log: %s\n", output);

	return fs;
}

void initShaders()
{
	// Create the programs

    gProgram[0] = glCreateProgram();
	gProgram[1] = glCreateProgram();

	// Create the shaders for both programs

    GLuint vs1 = createVS("vert.glsl");
    GLuint fs1 = createFS("frag.glsl");

	GLuint vs2 = createVS("vert2.glsl");
	GLuint fs2 = createFS("frag2.glsl");

	// Attach the shaders to the programs

	glAttachShader(gProgram[0], vs1);
	glAttachShader(gProgram[0], fs1);

	glAttachShader(gProgram[1], vs2);
	glAttachShader(gProgram[1], fs2);

	// Link the programs

    glLinkProgram(gProgram[0]);
	GLint status;
	glGetProgramiv(gProgram[0], GL_LINK_STATUS, &status);

	if (status != GL_TRUE)
	{
		cout << "Program link failed" << endl;
		exit(-1);
	}

	glLinkProgram(gProgram[1]);
	glGetProgramiv(gProgram[1], GL_LINK_STATUS, &status);

	if (status != GL_TRUE)
	{
		cout << "Program link failed" << endl;
		exit(-1);
	}

	// Get the locations of the uniform variables from both programs

	for (int i = 0; i < 2; ++i)
	{
		modelingMatrixLoc[i] = glGetUniformLocation(gProgram[i], "modelingMatrix");
		viewingMatrixLoc[i] = glGetUniformLocation(gProgram[i], "viewingMatrix");
		projectionMatrixLoc[i] = glGetUniformLocation(gProgram[i], "projectionMatrix");
		eyePosLoc[i] = glGetUniformLocation(gProgram[i], "eyePos");
	}
}


void initFlagVBO()
{
    GLuint vao;
    glGenVertexArrays(1, &vao);
    assert(vao > 0);
    glBindVertexArray(vao);
    cout << "vao = " << vao << endl;

	glEnableVertexAttribArray(0);
	glEnableVertexAttribArray(1);
	glEnableVertexAttribArray(2);
	assert(glGetError() == GL_NONE);

	glGenBuffers(1, &gVertexAttribBuffer);
	glGenBuffers(1, &gIndexBuffer);

	assert(gVertexAttribBuffer > 0 && gIndexBuffer > 0);

	glBindBuffer(GL_ARRAY_BUFFER, gVertexAttribBuffer);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gIndexBuffer);
	assert(glGetError() == GL_NONE);

///////////


	//glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
	//glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(gVertexDataSizeInBytes));
}

void initTexture()
{
    glGenTextures(1, &flagTexture);
    glBindTexture(GL_TEXTURE_2D, flagTexture);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_MIRRORED_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_MIRRORED_REPEAT);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    int width, height, nrChannels;
    rawImage = load_image(imagePath, &width, &height, &nrChannels);
    flagAspecRat = (float)width/height;

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, rawImage);
    glGenerateMipmap(GL_TEXTURE_2D);

    glUseProgram(gProgram[1]);
    glUniform1i(glGetUniformLocation(gProgram[1], "flagTexture"), 0);
    //std::cout << "get error: " << glGetError() << std::endl;

    free(rawImage);
}


void initFlag()
{
    //sampleBezier();

    glEnable(GL_DEPTH_TEST);
    initShaders();
    initFlagVBO();
    initTexture();
}



void drawBezier()
{
    Vertex controlPoints[4][4];
    gNormals = {} ;
    gVertices = {};
    gFaces = {};
    gTextures = {};

    static float wavingParam = 0;
    
    //controlPoints[0][0] = Vertex(0, 0,0);
    //controlPoints[0][1] = Vertex(0,-1,0);
    //controlPoints[0][2] = Vertex(0,-2,0);
    //controlPoints[0][3] = Vertex(0,-3,0);

    //controlPoints[1][0] = Vertex(3,0,1);
    //controlPoints[1][1] = Vertex(3,-1,2);
    //controlPoints[1][2] = Vertex(3,-2,2);
    //controlPoints[1][3] = Vertex(3,-3,1);

    //controlPoints[2][0] = Vertex(6,0,-1);
    //controlPoints[2][1] = Vertex(6,-1,-2);
    //controlPoints[2][2] = Vertex(6,-2,-2);
    //controlPoints[2][3] = Vertex(6,-3,-1);
    //
    //controlPoints[3][0] = Vertex(8,0,0);
    //controlPoints[3][1] = Vertex(8,-1,0);
    //controlPoints[3][2] = Vertex(8,-2,0);
    //controlPoints[3][3] = Vertex(8,-3,0);


    controlPoints[0][0] = Vertex(0, 0,0);
    controlPoints[0][1] = Vertex(0,-1,0);
    controlPoints[0][2] = Vertex(0,-2,0);
    controlPoints[0][3] = Vertex(0,-3,0);

    controlPoints[1][0] = Vertex(3,0,3 - wavingParam);
    controlPoints[1][1] = Vertex(3,-1,5 - wavingParam);
    controlPoints[1][2] = Vertex(3,-2,6 - wavingParam);
    controlPoints[1][3] = Vertex(3,-3,4 - wavingParam);

    controlPoints[2][0] = Vertex(6,0,-4 + wavingParam);
    controlPoints[2][1] = Vertex(6,-1,-6 + wavingParam);
    controlPoints[2][2] = Vertex(6,-2,-5 + wavingParam);
    controlPoints[2][3] = Vertex(6,-3,-3 + wavingParam);
    
    controlPoints[3][0] = Vertex(8,0,0);
    controlPoints[3][1] = Vertex(8,-1,0);
    controlPoints[3][2] = Vertex(8,-2,0);
    controlPoints[3][3] = Vertex(8,-3,0);

    static float changeWave = 0.15;
    if(wavingParam > 6 || wavingParam < 0)
    {
        changeWave *= -1;
    }
    wavingParam += changeWave;

    Vertex **Q; 

    Q = new Vertex*[sampleSize];
    GLfloat* samplePoints;

    samplePoints = linspace(0, 1, sampleSize);
    for(int row = 0; row < sampleSize; row++)
    {
        Q[row] = new Vertex[sampleSize];
        //Qy[row] = new Vertex[sampleSize];
        //Qz[row] = new Vertex[sampleSize];

        memset(&Q[row][0], 0, sizeof(struct Vertex)*sampleSize);
        //TODO remove logging debug
        //for(int col = 0; col< sampleSize; col++)
        //{
        //std::cout << "Q " << row << " " << col<<" init: " <<" "<< Q[row][col].x << " "<< Q[row][col].y <<" "<< Q[row][col].z << std::endl;
        //}
    }

    for(int s = 0; s < sampleSize; s++){
        for(int t = 0; t < sampleSize; t++){
            
            for(int i=0; i<=3; i++){
                for(int j = 0; j<=3; j++){
                    Q[s][t].x += bern(i, 3, samplePoints[s]) * bern(j, 3, samplePoints[t]) * controlPoints[i][j].x;
                    Q[s][t].y += bern(i, 3, samplePoints[s]) * bern(j, 3, samplePoints[t]) * controlPoints[i][j].y;
                    Q[s][t].z += bern(i, 3, samplePoints[s]) * bern(j, 3, samplePoints[t]) * controlPoints[i][j].z;
                    //TODO optimize bern function is called 3 times
                    //bern function can be computed more effiecntly
                }
            }

        }
    }

    
    for(int s = 0; s < sampleSize; s++){
        for(int t = 0; t < sampleSize; t++){
            gVertices.push_back(Q[s][t]);
        }
    }

    vector<glm::vec3> faceNormals; 
   
    //compute faces
    for(int s = 0; s <= sampleSize-2; s++)
    {
        for(int t = 0; t <= sampleSize-2; t++)
        {
            glm::vec3 v1, v2, v3, v4; 
            glm::vec3 norm1, norm2;

            int vIndex[3], tIndex[3];

            vIndex[0] = s*sampleSize + t;
            vIndex[1] = s*sampleSize + t + 1;
            vIndex[2] = (s+1)*sampleSize + t;

            //v1 = glm::vec3(Q[s][t].x, Q[s][t].y, Q[s][t].z);
            //v2 = glm::vec3(Q[s][t+1].x, Q[s][t+1].y, Q[s][t+1].z);
            //v3 = glm::vec3(Q[s+1][t].x, Q[s+1][t].y, Q[s+1][t].z);
            //v4 = glm::vec3(Q[s+1][t+1].x, Q[s+1][t+1].y, Q[s+1][t+1].z);

            v1 = Q[s][t].toGlmVec3();
            v2 = Q[s][t+1].toGlmVec3();
            v3 = Q[s+1][t].toGlmVec3();
            v4 = Q[s+1][t+1].toGlmVec3();

            //norm1 = cross(v2-v1, v3-v2);//TODO optimize
            //norm2 = cross(v3-v2, v4-v2);//v3-v2 is computed twice
            
            //faceNormals.push_back(norm1);
            //faceNormals.push_back(norm2);
            faceNormals.push_back(glm::triangleNormal(v1,v2,v3));
            faceNormals.push_back(glm::triangleNormal(v3,v2,v4));

            int ind = faceNormals.size()-2;
            glm::vec3 tempv = faceNormals[ind];
            //std::cout << "face normal " << ind << "" << tempv.x << " " << tempv.y <<" "<< tempv.z << std::endl;
            ind++;
            tempv = faceNormals[ind];
            //std::cout << "face normal " << ind << "" << tempv.x << " " << tempv.y <<" "<< tempv.z << std::endl;

            gFaces.push_back(Face(vIndex, vIndex, vIndex));//TODO critical
            vIndex[0] = (s+1)*sampleSize + t;
            vIndex[1] = s*sampleSize + t + 1;
            vIndex[2] = (s+1)*sampleSize + t + 1;

            gFaces.push_back(Face(vIndex, vIndex, vIndex));//TODO critical
        }
    }
    //faces ends


    gNormals.resize((sampleSize)*(sampleSize));

    gNormals[0] = Normal(faceNormals[0]);//top-left
    //std::cout << "normal " << faceNormals[(sampleSize-1)*2 -1].x << " "  << faceNormals[(sampleSize-1)*2 -1].y << " " << faceNormals[(sampleSize-1)*2 -1].z << " "<< std::endl;
    //std::cout << "normal " << faceNormals[(sampleSize-1 )*2 -2].x << " " << faceNormals[(sampleSize-1 )*2 -2].y << " " << faceNormals[(sampleSize-1 )*2 -2].z << " "<< std::endl;
    gNormals[sampleSize-1] = Normal((faceNormals[(sampleSize-1)*2 -1] + faceNormals[(sampleSize-1)*2 -2])/2.f);//bottom-left corner vertex

    gNormals[(sampleSize)*(sampleSize-1)] = Normal(
            (faceNormals[2*(sampleSize-2)*(sampleSize-1)] 
             + faceNormals[2*(sampleSize-2)*(sampleSize-1) + 1])/2.f
            );//top right

    gNormals[sampleSize*sampleSize -1] = Normal(
            faceNormals[2*(sampleSize-1) * (sampleSize-1)-1]);//bottom right

    for(int i = 1; i <= sampleSize-2; i++){
        int temp;
        //left edge of the patch
        gNormals[i] = Normal( (faceNormals[i*2 - 2] +
                    faceNormals[i*2 -1] + 
                    faceNormals[i*2]) / 3.f);

        //right edge of the patch
        temp = (sampleSize-1)*2*(sampleSize-2) + 2*i;
        gNormals[sampleSize*(sampleSize-1) + i] = Normal(
                    (faceNormals[temp-1] +
                     faceNormals[temp] +
                     faceNormals[temp+1]
                    )/3.f
                );
        //top edge vertex normals
        gNormals[sampleSize*i] = Normal(
                    (faceNormals[(sampleSize-1)*(i-1)] +
                     faceNormals[(sampleSize-1)*(i-1) + 1] +
                     faceNormals[(sampleSize-1) * i]
                    )/3.f
                );

        //bottom edge vertex normals
        gNormals[sampleSize*(i+1) - 1] = Normal(
                    (faceNormals[(sampleSize-1)*2*i -1] +
                     faceNormals[(sampleSize-1)*2*(i+1) -1] +
                     faceNormals[(sampleSize-1)*2*(i+1) -2]
                    )/3.f
                );

    }

    for(int s = 1; s <= sampleSize-2; s++)
    {
        for(int t = 1; t <= sampleSize-2; t++)
        {
            int temp;
            glm::vec3 norm;

            temp = ((sampleSize-1)*s + t) *2;
            norm = faceNormals[temp - 1] +
                        faceNormals[temp] +
                        faceNormals[temp +1];
            norm /= 3;
            gNormals[sampleSize*s + t] = Normal(norm);
        }

    }

    GLfloat u, v, step;
    u = 0;
    //v = 0;
    step = 1.0f/(sampleSize-1);
    for(int s = 0; s < sampleSize; s++)
    {//TODO HERE
        Texture tex;
        tex.u = u;
        v = 0;
        for(int t = 0; t < sampleSize; t++)
        {
            tex.v = v;
            gTextures.push_back(tex);
            v += step;
        }
        u += step;
    }

    



    //for(int patchNo = 0; patchNo < 1; patchNo++) 
    //{
    //    //gVertices.push_back(Vertex());

    //}


    delete[] samplePoints;
    for(int row = 0; row < sampleSize; row++)
    {
        delete[] Q[row];
    }
    
	gVertexDataSizeInBytes = gVertices.size() * 3 * sizeof(GLfloat);
	gNormalDataSizeInBytes = gNormals.size() * 3 * sizeof(GLfloat);
    gTextureDataSizeInBytes = gTextures.size() * 2 * sizeof(GLfloat);
	int indexDataSizeInBytes = gFaces.size() * 3 * sizeof(GLuint);
	GLfloat* vertexData = new GLfloat [gVertices.size() * 3];
	GLfloat* normalData = new GLfloat [gNormals.size() * 3];
	GLfloat* textureData = new GLfloat [gTextures.size() * 2];
	GLuint* indexData = new GLuint [gFaces.size() * 3];

    float minX = 1e6, maxX = -1e6;
    float minY = 1e6, maxY = -1e6;
    float minZ = 1e6, maxZ = -1e6;

    static int print_tex_flag = 1;
	for (int i = 0; i < gVertices.size(); ++i)
	{
        if(print_tex_flag){
            //std::cout << "vertex " << i << ": "<< gVertices[i].x << " "<< gVertices[i].y <<" "<< gVertices[i].z << std::endl;
        }

		vertexData[3*i] = gVertices[i].x;
		vertexData[3*i+1] = gVertices[i].y;
		vertexData[3*i+2] = gVertices[i].z;

        minX = std::min(minX, gVertices[i].x);
        maxX = std::max(maxX, gVertices[i].x);
        minY = std::min(minY, gVertices[i].y);
        maxY = std::max(maxY, gVertices[i].y);
        minZ = std::min(minZ, gVertices[i].z);
        maxZ = std::max(maxZ, gVertices[i].z);
	}

    //std::cout << "minX = " << minX << std::endl;
    //std::cout << "maxX = " << maxX << std::endl;
    //std::cout << "minY = " << minY << std::endl;
    //std::cout << "maxY = " << maxY << std::endl;
    //std::cout << "minZ = " << minZ << std::endl;
    //std::cout << "maxZ = " << maxZ << std::endl;

	for (int i = 0; i < gNormals.size(); ++i)
	{
        //std::cout << "normal " << i << " "<< gNormals[i].x << " "<< gNormals[i].y <<" "<< gNormals[i].z << std::endl;
		normalData[3*i] = gNormals[i].x;
		normalData[3*i+1] = gNormals[i].y;
		normalData[3*i+2] = gNormals[i].z;
	}

	for (int i = 0; i < gFaces.size(); ++i)
	{
		indexData[3*i] = gFaces[i].vIndex[0];
		indexData[3*i+1] = gFaces[i].vIndex[1];
		indexData[3*i+2] = gFaces[i].vIndex[2];
	}

    if(print_tex_flag){
        std::cout << "gTextures.size = " << gTextures.size() << std::endl;
    }
	for (int i = 0; i < gTextures.size(); ++i)//TODO optimize size calls
	{
		textureData[2*i] = gTextures[i].u;
		textureData[2*i+1] = gTextures[i].v;
        if(print_tex_flag){
            //std::cout << "tex << "<<i<<":  " << gTextures[i].u << " " << gTextures[i].v << std::endl;
        }
	}
	//for (int i = 0; i < gTextures.size(); ++i)//TODO optimize size calls
	//{
	//	textureData[2*i] = 0.f;
	//	textureData[2*i+1] = (float)i/gTextures.size();
    //    if(print_tex_flag){
    //        std::cout << "tex << "<<i<<":  \t" << textureData[2*i] << " " << textureData[2*i+1] << std::endl;
    //    }
    //}


	//glBufferData(GL_ARRAY_BUFFER, gVertexDataSizeInBytes + gNormalDataSizeInBytes, 0, GL_STATIC_DRAW);
	glBufferData(GL_ARRAY_BUFFER, gVertexDataSizeInBytes + gNormalDataSizeInBytes + gTextureDataSizeInBytes, 0, GL_STATIC_DRAW);
	glBufferSubData(GL_ARRAY_BUFFER, 0, gVertexDataSizeInBytes, vertexData);
	glBufferSubData(GL_ARRAY_BUFFER, gVertexDataSizeInBytes, gNormalDataSizeInBytes, normalData);
	glBufferSubData(GL_ARRAY_BUFFER, gVertexDataSizeInBytes + gNormalDataSizeInBytes, gTextureDataSizeInBytes, textureData);
	//glBufferSubData(GL_ARRAY_BUFFER, gVertexDataSizeInBytes + gNormalDataSizeInBytes, gTextureDataSizeInBytes, vertexData);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, indexDataSizeInBytes, indexData, GL_STATIC_DRAW);

    if(print_tex_flag){
        std::cout <<  "VertexDataSizeInBytes " <<  gVertexDataSizeInBytes << std::endl;
        std::cout <<  "gNormalDataSizeInBytes " <<  gNormalDataSizeInBytes << std::endl;
        std::cout <<  "gTextureDataSizeInBytes " <<  gTextureDataSizeInBytes << std::endl;
        std::cout <<  "get error " <<  glGetError() << std::endl;
    }

    print_tex_flag = 0;
	// done copying; can free now
	delete[] vertexData;
	delete[] normalData;
	delete[] indexData;
	delete[] textureData;
    
    glActiveTexture(GL_TEXTURE0 + 0);
    glBindTexture(GL_TEXTURE_2D, flagTexture);

	glBindBuffer(GL_ARRAY_BUFFER, gVertexAttribBuffer);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gIndexBuffer);

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(gVertexDataSizeInBytes));
	glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET((gVertexDataSizeInBytes + gNormalDataSizeInBytes)));


	glDrawElements(GL_TRIANGLES, gFaces.size() * 3, GL_UNSIGNED_INT, 0);
    static int should_draw = 1;
    if(should_draw){
        std::cout <<  "get error: " <<  glGetError() << std::endl;
    }
    should_draw = 0;
}

void genPatches()
{
    //for just one patch is generated
    gPatches = {};
    BezierPatch patch;

    static float wavingParam = 0;
    //patch.controlPoints[0][0] = Vertex(0, 0,0);
    //patch.controlPoints[0][1] = Vertex(0,-1,0);
    //patch.controlPoints[0][2] = Vertex(0,-2,0);
    //patch.controlPoints[0][3] = Vertex(0,-3,0);

    //patch.controlPoints[1][0] = Vertex(3,0,3 - wavingParam);
    //patch.controlPoints[1][1] = Vertex(3,-1,5 - wavingParam);
    //patch.controlPoints[1][2] = Vertex(3,-2,6 - wavingParam);
    //patch.controlPoints[1][3] = Vertex(3,-3,4 - wavingParam);

    //patch.controlPoints[2][0] = Vertex(6,0,-4 + wavingParam);
    //patch.controlPoints[2][1] = Vertex(6,-1,-6 + wavingParam);
    //patch.controlPoints[2][2] = Vertex(6,-2,-5 + wavingParam);
    //patch.controlPoints[2][3] = Vertex(6,-3,-3 + wavingParam);
    //
    //patch.controlPoints[3][0] = Vertex(8,0,0);
    //patch.controlPoints[3][1] = Vertex(8,-1,0);
    //patch.controlPoints[3][2] = Vertex(8,-2,0);
    //patch.controlPoints[3][3] = Vertex(8,-3,0);

    //patch.texTopLeft = Texture(0.f, 0.f);
    //patch.texTopRight = Texture(1.f, 0.f);
    //patch.texBottomLeft= Texture(0.f, 1.f);
    //patch.texBottomRight = Texture(1.f, 1.f);

    patch.controlPoints[0][0] = Vertex(0, 0,0);
    patch.controlPoints[0][1] = Vertex(0,-2,0);
    patch.controlPoints[0][2] = Vertex(0,-4,0);
    patch.controlPoints[0][3] = Vertex(0,-6,0);

    patch.controlPoints[1][0] = Vertex(1.25,0,3 - wavingParam);
    patch.controlPoints[1][1] = Vertex(1.25,-2,3 - wavingParam);
    patch.controlPoints[1][2] = Vertex(1.25,-4,3 - wavingParam);
    patch.controlPoints[1][3] = Vertex(1.25,-6,3 - wavingParam);

    patch.controlPoints[2][0] = Vertex(2.75,0,-3 + wavingParam);
    patch.controlPoints[2][1] = Vertex(2.75,-2,-3 + wavingParam);
    patch.controlPoints[2][2] = Vertex(2.75,-4,-3 + wavingParam);
    patch.controlPoints[2][3] = Vertex(2.75,-6,-3 + wavingParam);
    
    patch.controlPoints[3][0] = Vertex(4,0,0);
    patch.controlPoints[3][1] = Vertex(4,-2,0);
    patch.controlPoints[3][2] = Vertex(4,-4,0);
    patch.controlPoints[3][3] = Vertex(4,-6,0);

    patch.texTopLeft = Texture(0.f, 0.f);
    patch.texTopRight = Texture(0.5f, 0.f);
    patch.texBottomLeft= Texture(0.f, 1.0f);
    patch.texBottomRight = Texture(0.5f, 1.0f);

    gPatches.push_back(patch);

    patch.controlPoints[0][0] = Vertex(4, 0,0);
    patch.controlPoints[0][1] = Vertex(4,-2,0);
    patch.controlPoints[0][2] = Vertex(4,-4,0);
    patch.controlPoints[0][3] = Vertex(4,-6,0);

    patch.controlPoints[1][0] = patch.controlPoints[2][0].mirrorOver(patch.controlPoints[0][0]);
    patch.controlPoints[1][1] = patch.controlPoints[2][1].mirrorOver(patch.controlPoints[0][1]);
    patch.controlPoints[1][2] = patch.controlPoints[2][2].mirrorOver(patch.controlPoints[0][2]);
    patch.controlPoints[1][3] = patch.controlPoints[2][3].mirrorOver(patch.controlPoints[0][3]);


    patch.controlPoints[2][0] = Vertex(6.75,0,-3 + wavingParam);
    patch.controlPoints[2][1] = Vertex(6.75,-2,-3 + wavingParam);
    patch.controlPoints[2][2] = Vertex(6.75,-4,-3 + wavingParam);
    patch.controlPoints[2][3] = Vertex(6.75,-6,-3 + wavingParam);
    
    patch.controlPoints[3][0] = Vertex(8,0,0);
    patch.controlPoints[3][1] = Vertex(8,-2,0);
    patch.controlPoints[3][2] = Vertex(8,-4,0);
    patch.controlPoints[3][3] = Vertex(8,-6,0);

    patch.texTopLeft = Texture(0.5f, 0.f);
    patch.texTopRight = Texture(1.0f, 0.f);
    patch.texBottomLeft= Texture(0.5f, 1.0f);
    patch.texBottomRight = Texture(1.0f, 1.0f);

    gPatches.push_back(patch);

    static float changeWave = 0.15;
    if(wavingParam > 6 || wavingParam < 0)
    {
        changeWave *= -1;
    }
    wavingParam += changeWave;
}

int MAXwavyness =0;
void genSineWave()
{

    gPatches = {};
    int n, numOfPatches;
    n = patchSize;
    numOfPatches = n*n;

    static float time = 0;
    static float wavingParam = 0;
    static float changeWave = 0.025;
    static float pBound = 0.90f/log(2 + n);
    static float sin_2PI3 = sinf(M_PI*2.f/3.f);

    static glm::vec3 tl, bl, tr, br;//corners
    static GLfloat flagWidth = 8;
    static GLfloat flagHeight = flagWidth/flagAspecRat;
    flagWidth = 8;
    flagWidth *= 1.0f-(0.01*abs(wavingParam)/pBound);

    
    tl = glm::vec3(0.f, 0.f, 0.f);
    bl = glm::vec3(0.f, -flagHeight, 0.f);
    
    tr = glm::vec3(flagWidth, -flagHeight, 0.f);
    br = glm::vec3(flagWidth, -flagHeight, 0.f);


    GLfloat patchWidth = flagWidth / n;
    GLfloat patchHeight = flagHeight / n;

    GLfloat texSegLen = 1.0f/n;

    for(int i=0; i<n; ++i)
    {
        for(int j=0; j<n; ++j) {
            BezierPatch patch;
            GLfloat patchOffsetY = patchHeight*i;
            GLfloat patchOffsetX = patchWidth*j;

            patch.controlPoints[0][0] = Vertex(patchOffsetX, -(patchOffsetY ), 0);
            patch.controlPoints[0][1] = Vertex(patchOffsetX, -(patchOffsetY + patchHeight/3), 0);
            patch.controlPoints[0][2] = Vertex(patchOffsetX, -(patchOffsetY + 2 *patchHeight/3), 0);
            patch.controlPoints[0][3] = Vertex(patchOffsetX, -(patchOffsetY + patchHeight), 0);



            //TODO
            patchOffsetX += patchWidth/3;
            patch.controlPoints[1][0] = Vertex(patchOffsetX, -(patchOffsetY + sin_2PI3*wavingParam), 0+ sin_2PI3*wavingParam);
            patch.controlPoints[1][1] = Vertex(patchOffsetX, -(patchOffsetY + patchHeight/3+ sin_2PI3*wavingParam), 0+ sin_2PI3*wavingParam);
            patch.controlPoints[1][2] = Vertex(patchOffsetX, -(patchOffsetY + 2 *patchHeight/3+ sin_2PI3*wavingParam), 0+ sin_2PI3*wavingParam);
            patch.controlPoints[1][3] = Vertex(patchOffsetX, -(patchOffsetY + patchHeight+ sin_2PI3*wavingParam), 0+ sin_2PI3*wavingParam);

            //patch.controlPoints[1][0] = Vertex(patchOffsetX, -(patchOffsetY ), 0+ sin_2PI3*wavingParam);
            //patch.controlPoints[1][1] = Vertex(patchOffsetX, -(patchOffsetY + patchHeight/3), 0+ sin_2PI3*wavingParam);
            //patch.controlPoints[1][2] = Vertex(patchOffsetX, -(patchOffsetY + 2 *patchHeight/3), 0+ sin_2PI3*wavingParam);
            //patch.controlPoints[1][3] = Vertex(patchOffsetX, -(patchOffsetY + patchHeight), 0+ sin_2PI3*wavingParam);



            patchOffsetX += patchWidth/3;
            patch.controlPoints[2][0] = Vertex(patchOffsetX, -(patchOffsetY             - sin_2PI3*wavingParam), 0- sin_2PI3*wavingParam);
            patch.controlPoints[2][1] = Vertex(patchOffsetX, -(patchOffsetY + patchHeight/3 - sin_2PI3*wavingParam), 0- sin_2PI3*wavingParam);
            patch.controlPoints[2][2] = Vertex(patchOffsetX, -(patchOffsetY + 2 *patchHeight/3 - sin_2PI3*wavingParam), 0- sin_2PI3*wavingParam);
            patch.controlPoints[2][3] = Vertex(patchOffsetX, -(patchOffsetY + patchHeight - sin_2PI3*wavingParam), 0- sin_2PI3*wavingParam);
            
            //patch.controlPoints[2][0] = Vertex(patchOffsetX, -(patchOffsetY             ), 0- sin_2PI3*wavingParam);
            //patch.controlPoints[2][1] = Vertex(patchOffsetX, -(patchOffsetY + patchHeight/3 ), 0- sin_2PI3*wavingParam);
            //patch.controlPoints[2][2] = Vertex(patchOffsetX, -(patchOffsetY + 2 *patchHeight/3 ), 0- sin_2PI3*wavingParam);
            //patch.controlPoints[2][3] = Vertex(patchOffsetX, -(patchOffsetY + patchHeight ), 0- sin_2PI3*wavingParam);



            patchOffsetX += patchWidth/3;
            patch.controlPoints[3][0] = Vertex(patchOffsetX, -(patchOffsetY  ), 0);
            patch.controlPoints[3][1] = Vertex(patchOffsetX, -(patchOffsetY + patchHeight/3), 0);
            patch.controlPoints[3][2] = Vertex(patchOffsetX, -(patchOffsetY + 2 *patchHeight/3 ), 0);
            patch.controlPoints[3][3] = Vertex(patchOffsetX, -(patchOffsetY + patchHeight), 0);


            //for(int foo=0;)
            //{

            //}

            //////////// DO NOT TOUCH ()
            //patch.controlPoints[0][0] = Vertex(patchOffsetX, -(patchOffsetY + sin_2PI3*wavingParam), 0);
            //patch.controlPoints[0][1] = Vertex(patchOffsetX, -(patchOffsetY + patchHeight/3+ sin_2PI3*wavingParam), 0);
            //patch.controlPoints[0][2] = Vertex(patchOffsetX, -(patchOffsetY + 2 *patchHeight/3+ sin_2PI3*wavingParam), 0);
            //patch.controlPoints[0][3] = Vertex(patchOffsetX, -(patchOffsetY + patchHeight+ sin_2PI3*wavingParam), 0);



            ////TODO
            //patchOffsetX += patchWidth/3;
            //patch.controlPoints[1][0] = Vertex(patchOffsetX, -(patchOffsetY + sin_2PI3*wavingParam), 0);
            //patch.controlPoints[1][1] = Vertex(patchOffsetX, -(patchOffsetY + patchHeight/3+ sin_2PI3*wavingParam), 0);
            //patch.controlPoints[1][2] = Vertex(patchOffsetX, -(patchOffsetY + 2 *patchHeight/3+ sin_2PI3*wavingParam), 0);
            //patch.controlPoints[1][3] = Vertex(patchOffsetX, -(patchOffsetY + patchHeight+ sin_2PI3*wavingParam), 0);



            //patchOffsetX += patchWidth/3;
            //patch.controlPoints[2][0] = Vertex(patchOffsetX, -(patchOffsetY + sin_2PI3*wavingParam), 0);
            //patch.controlPoints[2][1] = Vertex(patchOffsetX, -(patchOffsetY + patchHeight/3+ sin_2PI3*wavingParam), 0);
            //patch.controlPoints[2][2] = Vertex(patchOffsetX, -(patchOffsetY + 2 *patchHeight/3+ sin_2PI3*wavingParam), 0);
            //patch.controlPoints[2][3] = Vertex(patchOffsetX, -(patchOffsetY + patchHeight+ sin_2PI3*wavingParam), 0);



            //patchOffsetX += patchWidth/3;
            //patch.controlPoints[3][0] = Vertex(patchOffsetX, -(patchOffsetY + sin_2PI3*wavingParam), 0);
            //patch.controlPoints[3][1] = Vertex(patchOffsetX, -(patchOffsetY + patchHeight/3+ sin_2PI3*wavingParam), 0);
            //patch.controlPoints[3][2] = Vertex(patchOffsetX, -(patchOffsetY + 2 *patchHeight/3+ sin_2PI3*wavingParam), 0);
            //patch.controlPoints[3][3] = Vertex(patchOffsetX, -(patchOffsetY + patchHeight+ sin_2PI3*wavingParam), 0);


            patch.texTopLeft = Texture(texSegLen*j, texSegLen*i);
            patch.texTopRight = Texture(texSegLen*(j+1), texSegLen*i);
            patch.texBottomLeft = Texture(texSegLen*j, texSegLen*(i+1));
            patch.texBottomRight = Texture(texSegLen*(j+1), texSegLen*(i+1));

            gPatches.push_back(patch);
        }
    }

    GLfloat sinVal = sin(time) * 0.25;
    time += 0.05;
    static int printtt =1;
    if(wavingParam > pBound || wavingParam < -pBound)
    {
        //if(printtt){
        //for(int i=0; i<gPatches.size(); i++){
        //    for(int foo=0; foo<4; ++foo){
        //        for(int bar=0; bar<4; ++bar){
        //            cout << "patch no: " <<i << " ctrlPT[" <<foo<< "][" << bar << "]    " <<
        //                "(" << gPatches[i].controlPoints[foo][bar].x << ", " <<  gPatches[i].controlPoints[foo][bar].y << ", " << gPatches[i].controlPoints[foo][bar].z << ") " << endl;
        //        }
        //    }
        //}
        //printtt = 0;
        //}
        //MAXwavyness = 1 && (printtt > 0);
        //printtt--;
        changeWave *= -1;
    }
    wavingParam += changeWave;
    //cout << "waving par: " << wavingParam << endl;

    //for(size_t i=0; i<)


}

void genSineWave2()
{
    gPatches = {}; 
    

    int dimCtrlPoints = patchSize*3 +1;
    Vertex** ctrlPoints = new Vertex*[dimCtrlPoints];
    for(int i=0; i < dimCtrlPoints; ++i)
    {
        ctrlPoints[i] = new Vertex[dimCtrlPoints];
    }

    for(int i=0; i<dimCtrlPoints; ++i)
    {
        for(int j=0; j<dimCtrlPoints; ++j)
        {
            //ctrlPoints[i][j] = 


        }
    }
    
        
    
    for(int i=0; i < dimCtrlPoints; ++i)
    {
        delete[] ctrlPoints[i];
    }
    delete[] ctrlPoints;
    
}

void drawFlag()
{
    gNormals = {} ;
    gVertices = {};
    gFaces = {};
    gTextures = {};

    //genPatches();
    genSineWave();
    
    //for(each patch sample patch)
    //;
    //draw


    GLfloat* samplePoints;
    samplePoints = linspace(0, 1, sampleSize);

    size_t numOfPatches = gPatches.size();
    for(size_t patchNo=0; patchNo < numOfPatches; ++patchNo)
    {
        //get control points
        vector<Vertex> patchVertices = {};
        vector<Normal> patchNormals = {} ;
        vector<Face> patchFaces = {};
        vector<Texture> patchTextures = {};

        vector<Vertex> fillQ(sampleSize, Vertex(0, 0, 0));
        vector<vector<Vertex>> Q(sampleSize, fillQ);//initialize with 0's

        for(int s = 0; s < sampleSize; s++)
        {
            for(int t = 0; t < sampleSize; t++)
            {

                for(int i=0; i<=3; i++)
                {
                    GLfloat bern_s = bern(i, 3, samplePoints[s]);
                    for(int j = 0; j<=3; j++)
                    {
                        GLfloat bern_t = bern(j, 3, samplePoints[t]);
                        Q[s][t].x += bern_s * bern_t * gPatches[patchNo].controlPoints[i][j].x;
                        Q[s][t].y += bern_s * bern_t * gPatches[patchNo].controlPoints[i][j].y;
                        Q[s][t].z += bern_s * bern_t * gPatches[patchNo].controlPoints[i][j].z;
                        //TODO optimize bern function is called 3 times
                        //bern function can be computed more effiecntly
                    }
                }

            }
        }

        //int temp_v=0;
        //static int should_print = patchSize;
        //if(patchNo%patchSize == 0 && MAXwavyness)
        //{
        //    for(int q=0; q<sampleSize; q++)
        //    {
        //        for(int r=0; r<sampleSize; r++)
        //        {
        //            cout << "vert "<< temp_v++ << ": "<< Q[q][r].x << ", "<< Q[q][r].y <<  ", "<< Q[q][r].z <<endl;
        //        }
        //    }
        //}

        for(int s = 0; s < sampleSize; s++)
        {
            for(int t = 0; t < sampleSize; t++)
            {
                patchVertices.push_back(Q[s][t]);
            }
        }


        //compute faces and their normals
        vector<glm::vec3> patchFaceNormals;

        size_t indexOffset = gVertices.size();

        for(int s = 0; s <= sampleSize-2; ++s)
        {
            for(int t = 0; t <= sampleSize-2; ++t)
            {
                glm::vec3 v1, v2, v3, v4;

                v1 = Q[s][t].toGlmVec3();
                v2 = Q[s][t+1].toGlmVec3();
                v3 = Q[s+1][t].toGlmVec3();
                v4 = Q[s+1][t+1].toGlmVec3();

                //compute normals
                patchFaceNormals.push_back(glm::triangleNormal(v1,v2,v3));
                patchFaceNormals.push_back(glm::triangleNormal(v3,v2,v4));

                //compute vertex indices for faces
                int vIndex[3];

                vIndex[0] = indexOffset + s*sampleSize + t;
                vIndex[1] = indexOffset + s*sampleSize + t + 1;
                vIndex[2] = indexOffset + (s+1)*sampleSize + t;
                patchFaces.push_back(Face(vIndex, vIndex, vIndex));

                vIndex[0] = indexOffset + (s+1)*sampleSize + t;
                //vIndex[1] = indexOffset + s*sampleSize + t + 1;//it is the same vertex
                vIndex[2] = indexOffset + (s+1)*sampleSize + t + 1;
                patchFaces.push_back(Face(vIndex, vIndex, vIndex));//TODO critical
            }
        }

        


        patchNormals.resize(sampleSize*sampleSize);

        //compute normals at corners
        patchNormals[0] = Normal(patchFaceNormals[0]);//top left
        patchNormals[sampleSize-1] = (patchFaceNormals[(sampleSize-1)*2 -1] + patchFaceNormals[(sampleSize-1)*2 -2])/2.f; //bottom-left
        patchNormals[(sampleSize)*(sampleSize-1)] = Normal(
                (patchFaceNormals[2*(sampleSize-2)*(sampleSize-1)] 
                 + patchFaceNormals[2*(sampleSize-2)*(sampleSize-1) + 1])/2.f
                );//top right


        patchNormals[sampleSize*sampleSize -1] = Normal(
                patchFaceNormals[2*(sampleSize-1) * (sampleSize-1)-1]);//bottom right

        //normals at the edges
        for(int i = 1; i <= sampleSize-2; i++)
        {
            int temp;
            //left edge of the patch
            patchNormals[i] = Normal( 
                    (patchFaceNormals[i*2 - 2] +
                     //patchFaceNormals[i*2 -1] + 
                     patchFaceNormals[i*2]
                    ) / 2.f
                    );

            //right edge of the patch
            temp = (sampleSize-1)*2*(sampleSize-2) + 2*i;
            patchNormals[sampleSize*(sampleSize-1) + i] = Normal(
                    (patchFaceNormals[temp-1] +
                     patchFaceNormals[temp]// +
                     //patchFaceNormals[temp+1]
                    )/2.f
                    );
            //top edge vertex normals
            patchNormals[sampleSize*i] = Normal(
                    (patchFaceNormals[(sampleSize-1)*2 *(i-1)] +
                     //patchFaceNormals[(sampleSize-1)*2 *(i-1) + 1] +
                     patchFaceNormals[(sampleSize-1)*2  * i]
                    )/2.f
                    );

            //bottom edge vertex normals
            patchNormals[sampleSize*(i+1) - 1] = Normal(
                    (patchFaceNormals[(sampleSize-1)*2*i -1] +
                     //patchFaceNormals[(sampleSize-1)*2*(i+1) -1] +
                     patchFaceNormals[(sampleSize-1)*2*(i+1) -2]
                    )/2.f
                    );

            //TODO fix
            ////left edge of the patch
            //patchNormals[i] = Normal( (patchFaceNormals[i*2 - 2] +
            //            patchFaceNormals[i*2 -1] + 
            //            patchFaceNormals[i*2]) / 3.f);

            ////right edge of the patch
            //temp = (sampleSize-1)*2*(sampleSize-2) + 2*i;
            //patchNormals[sampleSize*(sampleSize-1) + i] = Normal(
            //        (patchFaceNormals[temp-1] +
            //         patchFaceNormals[temp] +
            //         patchFaceNormals[temp+1]
            //        )/3.f
            //        );
            ////top edge vertex normals
            //patchNormals[sampleSize*i] = Normal(
            //        (patchFaceNormals[(sampleSize-1)*2 *(i-1)] +
            //         patchFaceNormals[(sampleSize-1)*2 *(i-1) + 1] +
            //         patchFaceNormals[(sampleSize-1)*2  * i]
            //        )/3.f
            //        );

            ////bottom edge vertex normals
            //patchNormals[sampleSize*(i+1) - 1] = Normal(
            //        (patchFaceNormals[(sampleSize-1)*2*i -1] +
            //         patchFaceNormals[(sampleSize-1)*2*(i+1) -1] +
            //         patchFaceNormals[(sampleSize-1)*2*(i+1) -2]
            //        )/3.f
            //        );

        }

        for(int s = 1; s <= sampleSize-2; s++)
        {
            for(int t = 1; t <= sampleSize-2; t++)
            {
                int temp;
                glm::vec3 norm;

                temp = ((sampleSize-1)*s + t) *2;
                //norm = patchFaceNormals[temp - 1] +
                //    patchFaceNormals[temp] +
                //    patchFaceNormals[temp +1];
                //norm /= 3;
                //TODO why are there only 3 normals in this computation
                //TODO FIX
                //TODO fix
                norm = patchFaceNormals[temp - 2] +
                        patchFaceNormals[temp -1] +
                        patchFaceNormals[temp];

                temp -= (sampleSize-1)*2;
                norm += patchFaceNormals[temp-1] +
                            patchFaceNormals[temp] +
                            patchFaceNormals[temp+1];
                norm /= 6.f;
                patchNormals[sampleSize*s + t] = Normal(norm);
            }

        }



        Texture topLeft, topRight, bottomLeft, bottomRight;
        topLeft = gPatches[patchNo].texTopLeft;  
        topRight = gPatches[patchNo].texTopRight;
        bottomLeft = gPatches[patchNo].texBottomLeft;
        bottomRight = gPatches[patchNo].texBottomRight;

        GLfloat u = 0.f;
        GLfloat v = 0.f;
        GLfloat step = 1.0f/(sampleSize-1);

        for(int s=0; s < sampleSize; s++)
        {
            Texture top = (topLeft*(1-u) + topRight * u );
            Texture bottom = (bottomLeft*(1-u) + bottomRight * u);
            
            v = 0.f;
            for(int t=0; t < sampleSize; ++t)
            {
                Texture tex = (top*(1-v) + bottom* v);

                patchTextures.push_back(tex);

                v += step;
            }
            u +=step;
        }
        
        size_t numOfVertices = patchVertices.size();
        for(size_t i=0; i<numOfVertices; ++i)
        {
            gVertices.push_back(patchVertices[i]);
        }

        size_t numOfNormals = patchNormals.size();
        for(size_t i=0; i<numOfNormals; ++i)
        {
            gNormals.push_back(patchNormals[i]);
        }

        size_t numOfFaces = patchFaces.size();
        for(size_t i=0; i<numOfFaces; ++i)
        {
            gFaces.push_back(patchFaces[i]);
        }

        size_t numOfTextures = patchTextures.size();
        for(size_t i=0; i< numOfTextures; ++i)
        {
            gTextures.push_back(patchTextures[i]);
        }

        //pacth loop ends
    }

    // normals of vertices at the patch edges are incorrect fix it

    //for(int patchNo; patchNo < numOfPatches)

    delete[] samplePoints;

	gVertexDataSizeInBytes = gVertices.size() * 3 * sizeof(GLfloat);
	gNormalDataSizeInBytes = gNormals.size() * 3 * sizeof(GLfloat);
    gTextureDataSizeInBytes = gTextures.size() * 2 * sizeof(GLfloat);
	int indexDataSizeInBytes = gFaces.size() * 3 * sizeof(GLuint);
	GLfloat* vertexData = new GLfloat [gVertices.size() * 3];
	GLfloat* normalData = new GLfloat [gNormals.size() * 3];
	GLfloat* textureData = new GLfloat [gTextures.size() * 2];
	GLuint* indexData = new GLuint [gFaces.size() * 3];

    float minX = 1e6, maxX = -1e6;
    float minY = 1e6, maxY = -1e6;
    float minZ = 1e6, maxZ = -1e6;

    static int print_tex_flag = 1;
	for (int i = 0; i < gVertices.size(); ++i)
	{
        if(print_tex_flag){
            //std::cout << "vertex " << i << ": "<< gVertices[i].x << " "<< gVertices[i].y <<" "<< gVertices[i].z << std::endl;
        }

		vertexData[3*i] = gVertices[i].x;
		vertexData[3*i+1] = gVertices[i].y;
		vertexData[3*i+2] = gVertices[i].z;

        minX = std::min(minX, gVertices[i].x);
        maxX = std::max(maxX, gVertices[i].x);
        minY = std::min(minY, gVertices[i].y);
        maxY = std::max(maxY, gVertices[i].y);
        minZ = std::min(minZ, gVertices[i].z);
        maxZ = std::max(maxZ, gVertices[i].z);
	}

    //std::cout << "minX = " << minX << std::endl;
    //std::cout << "maxX = " << maxX << std::endl;
    //std::cout << "minY = " << minY << std::endl;
    //std::cout << "maxY = " << maxY << std::endl;
    //std::cout << "minZ = " << minZ << std::endl;
    //std::cout << "maxZ = " << maxZ << std::endl;

	for (int i = 0; i < gNormals.size(); ++i)
	{
        //std::cout << "normal " << i << " "<< gNormals[i].x << " "<< gNormals[i].y <<" "<< gNormals[i].z << std::endl;
		normalData[3*i] = gNormals[i].x;
		normalData[3*i+1] = gNormals[i].y;
		normalData[3*i+2] = gNormals[i].z;
	}

	for (int i = 0; i < gFaces.size(); ++i)
	{
		indexData[3*i] = gFaces[i].vIndex[0];
		indexData[3*i+1] = gFaces[i].vIndex[1];
		indexData[3*i+2] = gFaces[i].vIndex[2];
	}

    if(print_tex_flag){
        std::cout << "gTextures.size = " << gTextures.size() << std::endl;
    }
	for (int i = 0; i < gTextures.size(); ++i)//TODO optimize size calls
	{
		textureData[2*i] = gTextures[i].u;
		textureData[2*i+1] = gTextures[i].v;
        if(print_tex_flag){
            //std::cout << "tex << "<<i<<":  " << gTextures[i].u << " " << gTextures[i].v << std::endl;
        }
	}
	//for (int i = 0; i < gTextures.size(); ++i)//TODO optimize size calls
	//{
	//	textureData[2*i] = 0.f;
	//	textureData[2*i+1] = (float)i/gTextures.size();
    //    if(print_tex_flag){
    //        std::cout << "tex << "<<i<<":  \t" << textureData[2*i] << " " << textureData[2*i+1] << std::endl;
    //    }
    //}


	//glBufferData(GL_ARRAY_BUFFER, gVertexDataSizeInBytes + gNormalDataSizeInBytes, 0, GL_STATIC_DRAW);
	glBufferData(GL_ARRAY_BUFFER, gVertexDataSizeInBytes + gNormalDataSizeInBytes + gTextureDataSizeInBytes, 0, GL_STATIC_DRAW);
	glBufferSubData(GL_ARRAY_BUFFER, 0, gVertexDataSizeInBytes, vertexData);
	glBufferSubData(GL_ARRAY_BUFFER, gVertexDataSizeInBytes, gNormalDataSizeInBytes, normalData);
	glBufferSubData(GL_ARRAY_BUFFER, gVertexDataSizeInBytes + gNormalDataSizeInBytes, gTextureDataSizeInBytes, textureData);
	//glBufferSubData(GL_ARRAY_BUFFER, gVertexDataSizeInBytes + gNormalDataSizeInBytes, gTextureDataSizeInBytes, vertexData);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, indexDataSizeInBytes, indexData, GL_STATIC_DRAW);

    if(print_tex_flag){
        std::cout <<  "VertexDataSizeInBytes " <<  gVertexDataSizeInBytes << std::endl;
        std::cout <<  "gNormalDataSizeInBytes " <<  gNormalDataSizeInBytes << std::endl;
        std::cout <<  "gTextureDataSizeInBytes " <<  gTextureDataSizeInBytes << std::endl;
        std::cout <<  "get error " <<  glGetError() << std::endl;
    }

    print_tex_flag = 0;
	// done copying; can free now
	delete[] vertexData;
	delete[] normalData;
	delete[] indexData;
	delete[] textureData;
    
    glActiveTexture(GL_TEXTURE0 + 0);
    glBindTexture(GL_TEXTURE_2D, flagTexture);

	glBindBuffer(GL_ARRAY_BUFFER, gVertexAttribBuffer);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gIndexBuffer);

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(gVertexDataSizeInBytes));
	glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET((gVertexDataSizeInBytes + gNormalDataSizeInBytes)));


	glDrawElements(GL_TRIANGLES, gFaces.size() * 3, GL_UNSIGNED_INT, 0);
    static int should_draw = 1;
    if(should_draw){
        std::cout <<  "get error: " <<  glGetError() << std::endl;
    }
    should_draw = 0;

}//draw flag ends

void display()
{
    glClearColor(0, 0, 0, 1);
    glClearDepth(1.0f);
    glClearStencil(0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

	static float angle = 0;

	float angleRad = (float) (angle / 180.0) * M_PI;
	
	// Compute the modeling matrix

	////modelingMatrix = glm::translate(glm::mat4(1.0), glm::vec3(0.0f, -0.4f, -5.0f));
	////modelingMatrix = glm::rotate(modelingMatrix, angleRad, glm::vec3(0.0, 1.0, 0.0));
    //glm::mat4 matT = glm::translate(glm::mat4(1.0), glm::vec3(-0.5f, -0.4f, -5.0f));   // same as above but more clear//original AOA's version TODO
    ////glm::mat4 matR = glm::rotate(glm::mat4(1.0), angleRad, glm::vec3(0.0, 1.0, 0.0));
    //glm::mat4 matRx = glm::rotate<float>(glm::mat4(1.0), (-90. / 180.) * M_PI, glm::vec3(1.0, 0.0, 0.0));
    //glm::mat4 matRy = glm::rotate<float>(glm::mat4(1.0), (-90. / 180.) * M_PI, glm::vec3(0.0, 1.0, 0.0));
    //glm::mat4 matRz = glm::rotate<float>(glm::mat4(1.0), angleRad, glm::vec3(0.0, 0.0, 1.0));

    //modelingMatrix = matRy * matRx;

    static float wavingAngle = 0.f;

    glm::mat4 matRx = glm::rotate<float>(glm::mat4(1.0), (-12. / 180.) * M_PI, glm::vec3(1.0, 0.0, 0.0));
    glm::mat4 matRy = glm::rotate<float>(glm::mat4(1.0), (wavingAngle / 180.) * M_PI, glm::vec3(0.0, 1.0, 0.0));
    //glm::mat4 matRz = glm::rotate<float>(glm::mat4(1.0), angleRad, glm::vec3(0.0, 0.0, 1.0));
    glm::mat4 matT = glm::translate(glm::mat4(1.0), glm::vec3(-4.f, 1.5f, -15.0f));   // same as above but more clear//original AOA's version TODO
    modelingMatrix = matT * matRx * matRy;

    const static float wavingBound = 5;
    static float changeWave = 0.1;
    wavingAngle += changeWave;
    if(wavingAngle > wavingBound || wavingAngle < -wavingBound)
    {
        changeWave *= -1;
    }

    // Let's make some alternating roll rotation
    static float rollDeg = 0;
    //static float changeRoll = 2.5;
    static float changeRoll = 0.25;
    float rollRad = (float) (rollDeg / 180.f) * M_PI;
    //rollDeg += changeRoll;
    if (rollDeg >= 10.f || rollDeg <= -10.f)
    {
        changeRoll *= -1.f;
    }
    glm::mat4 matRoll = glm::rotate<float>(glm::mat4(1.0), rollRad, glm::vec3(1.0, 0.0, 0.0));

    // Let's make some pitch rotation
    static float pitchDeg = 0;
    static float changePitch = 0.1;
    float startPitch = 0;
    float endPitch = 90;
    float pitchRad = (float) (pitchDeg / 180.f) * M_PI;
    //pitchDeg += changePitch;
    if (pitchDeg >= endPitch)
    {
        changePitch = 0;
    }
    //glm::mat4 matPitch = glm::rotate<float>(glm::mat4(1.0), pitchRad, glm::vec3(0.0, 0.0, 1.0));
    //modelingMatrix = matRoll * matPitch * modelingMatrix; // gimbal lock
    //modelingMatrix = matPitch * matRoll * modelingMatrix;   // no gimbal lock

    glm::quat q0(0, 1, 0, 0); // along x
    glm::quat q1(0, 0, 1, 0); // along y
    glm::quat q = glm::mix(q0, q1, (pitchDeg - startPitch) / (endPitch - startPitch));

    float sint = sin(rollRad / 2);
    glm::quat rollQuat(cos(rollRad/2), sint * q.x, sint * q.y, sint * q.z);
    glm::quat pitchQuat(cos(pitchRad/2), 0, 0, 1 * sin(pitchRad/2));
    //modelingMatrix = matT * glm::toMat4(pitchQuat) * glm::toMat4(rollQuat) * modelingMatrix;
    //modelingMatrix = matT * glm::toMat4(rollQuat) * glm::toMat4(pitchQuat) * modelingMatrix; // roll is based on pitch
    //modelingMatrix =  modelingMatrix; // roll is based on pitch

    //cout << rollQuat.w << " " << rollQuat.x << " " << rollQuat.y << " " << rollQuat.z << endl;

	// Set the active program and the values of its uniform variables

	glUseProgram(gProgram[activeProgramIndex]);
	glUniformMatrix4fv(projectionMatrixLoc[activeProgramIndex], 1, GL_FALSE, glm::value_ptr(projectionMatrix));
	glUniformMatrix4fv(viewingMatrixLoc[activeProgramIndex], 1, GL_FALSE, glm::value_ptr(viewingMatrix));
	glUniformMatrix4fv(modelingMatrixLoc[activeProgramIndex], 1, GL_FALSE, glm::value_ptr(modelingMatrix));
	glUniform3fv(eyePosLoc[activeProgramIndex], 1, glm::value_ptr(eyePos));

	// Draw the scene
    //drawModel();
    //drawBezier();
    drawFlag();

	//angle += 0.5;
}

void reshape(GLFWwindow* window, int w, int h)
{
    w = w < 1 ? 1 : w;
    h = h < 1 ? 1 : h;

    gWidth = w;
    gHeight = h;

    glViewport(0, 0, w, h);

    //glMatrixMode(GL_PROJECTION);
    //glLoadIdentity();
    //glOrtho(-10, 10, -10, 10, -10, 10);
    //gluPerspective(45, 1, 1, 100);

	// Use perspective projection

	float fovyRad = (float) (45.0 / 180.0) * M_PI;
	projectionMatrix = glm::perspective(fovyRad, 1.0f, 1.0f, 100.0f);

	// Assume default camera position and orientation (camera is at
	// (0, 0, 0) with looking at -z direction and its up vector pointing
	// at +y direction)

	viewingMatrix = glm::mat4(1);

    //glMatrixMode(GL_MODELVIEW);
    //glLoadIdentity();
}

void keyboard(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_Q && action == GLFW_PRESS)
    {
        glfwSetWindowShouldClose(window, GLFW_TRUE);
    }
    else if (key == GLFW_KEY_G && action == GLFW_PRESS)
    {
        //glShadeModel(GL_SMOOTH);
        activeProgramIndex = 0;
    }
    else if (key == GLFW_KEY_P && action == GLFW_PRESS)
    {
        //glShadeModel(GL_SMOOTH);
        activeProgramIndex = 1;
    }
    else if (key == GLFW_KEY_F && action == GLFW_PRESS)
    {
        //glShadeModel(GL_FLAT);
    }
    if(key == GLFW_KEY_W && action == GLFW_PRESS)
    {
        sampleSize++;
    }
    else if(key == GLFW_KEY_S && action == GLFW_PRESS && sampleSize > 2)
    {
        sampleSize--;
    }
    if(key == GLFW_KEY_E && action == GLFW_PRESS)
    {
        patchSize++;
    }
    else if(key == GLFW_KEY_D && action == GLFW_PRESS && patchSize > 1)
    {
        patchSize--;
    }


}

void mainLoop(GLFWwindow* window)
{
    while (!glfwWindowShouldClose(window))
    {
        display();
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
}

int main(int argc, char** argv)   // Create Main Function For Bringing It All Together
{
    GLFWwindow* window;

    if(argc > 1)
    {
        imagePath = strdup(argv[1]);
    }
    else
    {
        imagePath = strdup("metu_flag.jpg");
    }


    if (!glfwInit())
    {
        exit(-1);
    }

    //glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
    //glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, GL_TRUE);//TODO remove debug
    //glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_COMPAT_PROFILE);
    //glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    int width = 640, height = 480;
    window = glfwCreateWindow(width, height, "Simple Example", NULL, NULL);

    if (!window)
    {
        glfwTerminate();
        exit(-1);
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    // Initialize GLEW to setup the OpenGL Function pointers
    if (GLEW_OK != glewInit())
    {
        std::cout << "Failed to initialize GLEW" << std::endl;
        return EXIT_FAILURE;
    }

    char rendererInfo[512] = {0};
    strcpy(rendererInfo, (const char*) glGetString(GL_RENDERER));
    strcat(rendererInfo, " - ");
    strcat(rendererInfo, (const char*) glGetString(GL_VERSION));
    glfwSetWindowTitle(window, rendererInfo);

    //init();
    initFlag();

    glfwSetKeyCallback(window, keyboard);
    glfwSetWindowSizeCallback(window, reshape);

    reshape(window, width, height); // need to call this once ourselves
    mainLoop(window); // this does not return unless the window is closed

    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}
