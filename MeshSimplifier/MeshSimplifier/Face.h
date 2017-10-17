/*
Class defining a single face of the mesh
contains list of its vertices, its normal, color, and whether it is being used currently in the mesh (not used if mesh is simplified and face is destroyed)
*/
#pragma once
#include <iostream>
#include "Dependencies/glm/glm.hpp"
#include "Dependencies/glew/glew.h"
#include "Dependencies/glut/glut.h"
#include <vector>
#include "Vertex.h"
#include "Dependencies/glm/ext.hpp"

class Face
{

private:
	int index; //this face's index in the meshobject list
	std::vector<Vertex*> vertices;
	vec3 normal;
	vec3 Color;
	int vectorIndex;
	bool isUsed;

public:
	//error quadric of the face (used in quadric simplification)
	mat4x4 K;
	Face(Vertex*& one, Vertex*& two, Vertex*& three);
	//method to check if triangle is degenerate (has 2 same points)
	bool isDegenerate();
	//method to set a new vertex for the face when one vertex is deleted
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