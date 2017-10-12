#include "Vertex.h"
#include "Face.h"

bool Face::isDegenerate() {
	vec3 v0Position = (vertices[0])->position;
	vec3 v1Position = (vertices[1])->position;
	vec3 v2Position = (vertices[2])->position;
	if (v0Position[0] == v1Position[0] && v0Position[1] == v1Position[1] && v0Position[2] == v1Position[2]) return true;
	if (v0Position[0] == v2Position[0] && v0Position[1] == v2Position[1] && v0Position[2] == v2Position[2]) return true;
	if (v2Position[0] == v1Position[0] && v2Position[1] == v1Position[1] && v2Position[2] == v1Position[2]) return true;
	return false;
}

void Face::setNewVertForFace(Vertex* const oldVert, Vertex* const& newVert) {
	for (unsigned int i = 0; i < vertices.size(); i++) {
		if ((vertices[i])->position[0] == oldVert->position[0] &&
			(vertices[i])->position[1] == oldVert->position[1] &&
			(vertices[i])->position[2] == oldVert->position[2]) {
			vertices[i] = newVert;
			return;
		}
	}
	updateFace();
}

void Face::updateFace() {
	vec3 V = vertices[1]->position - vertices[0]->position;
	vec3 W = vertices[2]->position - vertices[0]->position;
	normal = normalize(glm::cross(V, W));
	float a = normal[0];
	float b = normal[1];
	float c = normal[2];
	float d = glm::dot(-(vertices[0]->position), normal);
	K[0] = { a * a, a * b, a * c, a * d };
	K[1] = { b * a, b * b, b * c, b * d };
	K[2] = { c * a, c * b, c * c, c * d };
	K[3] = { d * a, d * b, d * c, d * d };
}

void Face::setNotUsed() {
	isUsed = false;
	//vertices[0]->isUsed = false;
	//vertices[1]->isUsed = false;
	//vertices[2]->isUsed = false;
}
