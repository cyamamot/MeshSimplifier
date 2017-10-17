#include "Vertex.h"
#include "Face.h"

Face::Face(Vertex*& one, Vertex*& two, Vertex*& three)
{
	isUsed = true;
	vertices.push_back(one);
	vertices.push_back(two);
	vertices.push_back(three);
	vectorIndex = 0;
	updateFace();
}

bool Face::isDegenerate()
{
	vec3 v0Position = (vertices[0])->getPosition();
	vec3 v1Position = (vertices[1])->getPosition();
	vec3 v2Position = (vertices[2])->getPosition();
	//checks if two faces are in the exact same position and location
	if (v0Position[0] == v1Position[0] && v0Position[1] == v1Position[1] && v0Position[2] == v1Position[2]) return true;
	if (v0Position[0] == v2Position[0] && v0Position[1] == v2Position[1] && v0Position[2] == v2Position[2]) return true;
	if (v2Position[0] == v1Position[0] && v2Position[1] == v1Position[1] && v2Position[2] == v1Position[2]) return true;
	return false;
}

void Face::setNewVertForFace(Vertex* const oldVert, Vertex* const& newVert)
{
	for (unsigned int i = 0; i < vertices.size(); i++) 
	{
		if ((vertices[i])->getPosition()[0] == oldVert->getPosition()[0] &&
			(vertices[i])->getPosition()[1] == oldVert->getPosition()[1] &&
			(vertices[i])->getPosition()[2] == oldVert->getPosition()[2]) {
			vertices[i] = newVert;
			return;
		}
	}
	updateFace();
}

void Face::updateFace() 
{
	vec3 V = vertices[1]->getPosition() - vertices[0]->getPosition();
	vec3 W = vertices[2]->getPosition() - vertices[0]->getPosition();
	//recalculates normal for face if vertices change
	normal = normalize(glm::cross(V, W));
	float a = normal[0];
	float b = normal[1];
	float c = normal[2];
	float d = glm::dot(-(vertices[0]->getPosition()), normal);
	//racalculates this face's error quadric based on new comprising vertices
	K[0] = { a * a, a * b, a * c, a * d };
	K[1] = { b * a, b * b, b * c, b * d };
	K[2] = { c * a, c * b, c * c, c * d };
	K[3] = { d * a, d * b, d * c, d * d };
}

void Face::setNotUsed()
{
	isUsed = false;
}

int Face::getIndex()
{
	return index;
}

void Face::setIndex(int ind) 
{
	index = ind;
}

Vertex* Face::getVertex(int ind)
{
	return vertices[ind];
}

void Face::setVertex(int ind, Vertex* vert)
{
	vertices[ind] = vert;
}

vec3 Face::getNormal()
{
	return normal;
}

void Face::setNormal(vec3 norm)
{
	normal = norm;
}

vec3 Face::getColor() 
{
	return Color;
}

void Face::setColor(vec3 col)
{
	Color = col;
}

int Face::getVectorIndex()
{
	return vectorIndex;
}

void Face::setVectorIndex(int vecInd)
{
	vectorIndex = vecInd;
}

bool Face::getIsUsed() 
{
	return isUsed;
}

void Face::setIsUsed(bool use)
{
	isUsed = use;
}