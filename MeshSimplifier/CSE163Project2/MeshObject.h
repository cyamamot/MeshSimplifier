#pragma once

#include <iostream>
#include "Dependencies/glm/glm.hpp"
#include "Dependencies/glew/glew.h"
#include "Dependencies/glut/glut.h"
#include <vector>
#include "Vertex.h"
#include "Face.h"
#include <list>
#include <utility>

typedef glm::vec4 vec4;

class MeshObject {

private:
	int numVertices;
	int numFaces;
	std::vector<Vertex*> vertList; //list of pointers to vertices
	std::vector<Face*> faceList; //list of pointers to faces
	std::vector<vec3> vertexPositionList; //list of the positions of every vertex in order
	std::vector<vec3> normalsList; //list of normals of faces in order
	std::vector<vec3> faceColorList; //list of color of each face
	std::vector<int> vertexIndexList; //list specifying order of indices in prior lists 
	std::list<int> unsimplificationOrder;

	double Xmin = 100000.0, Ymin = 100000.0, Zmin = 100000.0;
	double Xmax = -100000.0, Ymax = -100000.0, Zmax = -100000.0;

	//number of times to collapse vertex pairs based on quadratic simplification
	int simplificationIterations = 25;

public:
	// Default associated variables
	GLuint modelviewPos;
	glm::mat4 modelview;
	//default constructor
	MeshObject(const char* filename);
	void setVertices(std::vector<float>& v, std::ifstream& infile);
	void setFaces(std::vector<float>& v, std::ifstream& infile);
	int getNumFaces();
	float norm();
	void initialize(GLuint &VAO, GLuint &VBO, GLuint &EBO, GLuint &NBO, GLuint &CBO);
	void drawObject(GLuint &VAO, GLuint &VBO, GLuint &EBO, GLuint &NBO, GLuint &CBO);
	void setVertexNormals();
	void edgeCollapse(int vert1, int vert2);
	void undoCollapse(int vert1, int vert2);
	void updateVBO(GLuint &VAO, GLuint &VBO, GLuint &EBO, GLuint &NBO, GLuint &CBO);
	struct Comp;
	struct quadricContainer;
	void quadricSimplification();
	void quadricUnsimplification();
	void printMatrix(mat4x4& m);
};


