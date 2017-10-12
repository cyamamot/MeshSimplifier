#pragma once
#include <iostream>
#include "Dependencies/glm/glm.hpp"
#include "Dependencies/glew/glew.h"
#include "Dependencies/glut/glut.h"
#include <vector>
#include "Vertex.h"
#include "Dependencies/glm/ext.hpp"

class Face {

private:
	int index; //this face's index in the meshobject list
	std::vector<Vertex*> vertices;
	vec3 normal;
	vec3 Color;
	int vectorIndex;
	bool isUsed;

public:
	mat4x4 K;
	Face(Vertex*& one, Vertex*& two, Vertex*& three);
	//method to check if triangle is degenerate (has 2 same points)
	bool isDegenerate();
	//method to set the old collapsed vertex as the newly creade midpoint one
	void setNewVertForFace(Vertex* const oldVert, Vertex* const& newVert);
	//update member variables
	void updateFace();
	void setNotUsed();
	//list of getters and setters
	int getIndex();
	void setIndex(int ind);
	Vertex* getVertex(int ind);
	void setVertex(int ind, Vertex* vert);
	vec3 getNormal();
	void setNormal(vec3 norm);
	vec3 getColor();
	void setColor(vec3 col);
	int getVectorIndex();
	void setVectorIndex(int vecInd);
	bool getIsUsed();
	void setIsUsed(bool use);
};