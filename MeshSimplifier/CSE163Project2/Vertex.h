#pragma once
#include <iostream>
#include "Dependencies/glm/glm.hpp"
#include "Dependencies/glew/glew.h"
#include "Dependencies/glut/glut.h"
#include <vector>
#include <algorithm>  

class Face;

typedef glm::vec3 vec3;
typedef glm::mat4x4 mat4x4;

class Vertex {

private:
	int index; //this vertex's index in the meshobject list
	vec3 position;
	vec3 normal;
	std::vector<Vertex*> adjacentVertices;
	std::vector<Face*> adjacentFaces;
	bool isUsed;
	int vertexVectorIndex;
	int faceVectorIndex;
	Vertex* replacedBy;

public:
	mat4x4 Q;

	Vertex(float x, float y, float z) {
		isUsed = true;
		position = vec3(x, y, z);
		vertexVectorIndex = 0;
		faceVectorIndex = 0;
	}
	//getters and setters
	int getIndex();
	void setIndex(int ind);
	vec3 getPosition();
	void setPosition(vec3 pos);
	vec3 getNormal();
	Face* getFace(int ind);
	int getAdjacentFaceListSize();
	Vertex* getVertex(int ind);
	int getAdjacentVertexListSize();
	void pushFace(Face* face);
	void pushVertex(Vertex* vert);
	bool getIsUsed();
	void setIsUsed(bool set);
	void setReplacedBy(Vertex* rep);
	//Let new vertex absorb adjacencies of old vertices
	void absorb(Vertex*& vert1, Vertex*& vert2);
	//method to set the old collapsed vertex as the newly create midpoint one
	void setNewVertForVert(Vertex* const oldVert, Vertex* const& newVert);
	//method that sets this in all adjacent vertices and faces as the newly created midpoint vertex
	void collapsed(Vertex*& mid);
	void revertVertex();
	void cleanMesh();
	void removeDegenerate();
	void updateVertex();
	void setNormal();
};