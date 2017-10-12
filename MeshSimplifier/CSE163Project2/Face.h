#pragma once
#include <iostream>
#include "Dependencies/glm/glm.hpp"
#include "Dependencies/glew/glew.h"
#include "Dependencies/glut/glut.h"
#include <vector>
#include "Vertex.h"
#include "Dependencies/glm/ext.hpp"

class Face {

public:
	int index; //this face's index in the meshobject list
	std::vector<Vertex*> vertices;
	vec3 normal;
	vec3 Color;
	int vectorIndex;
	bool isUsed;

	mat4x4 K;

	Face(Vertex*& one, Vertex*& two, Vertex*& three) {
		isUsed = true;
		vertices.push_back(one);
		vertices.push_back(two);
		vertices.push_back(three);
		vectorIndex = 0;
		updateFace();
	}

	//method to check if triangle is degenerate (has 2 same points)
	bool isDegenerate();

	//method to set the old collapsed vertex as the newly creade midpoint one
	void setNewVertForFace(Vertex* const oldVert, Vertex* const& newVert);

	//update member variables
	void updateFace();

	void setNotUsed();
};