#include "Vertex.h"
#include "Face.h"

void Vertex::absorb(Vertex*& vert1, Vertex*& vert2) {
	for (unsigned int i = 0; i < vert1->adjacentVertices.size(); i++) {
		Vertex* temp1 = vert1->adjacentVertices[i];
		if (temp1 == vert2 || temp1 == this) continue;
		adjacentVertices.push_back(vert1->adjacentVertices[i]);
	}
	for (unsigned int j = 0; j < vert2->adjacentVertices.size(); j++) {
		Vertex* temp2 = vert2->adjacentVertices[j];
		if (temp2 == vert1 || temp2 == this) continue;
		adjacentVertices.push_back(vert2->adjacentVertices[j]);
	}
	for (unsigned int k = 0; k < vert1->adjacentFaces.size(); k++) {
		adjacentFaces.push_back(vert1->adjacentFaces[k]);
	}
	for (unsigned int l = 0; l < vert2->adjacentFaces.size(); l++) {
		adjacentFaces.push_back(vert2->adjacentFaces[l]);
	}
}

void Vertex::setNewVertForVert(Vertex* const oldVert, Vertex* const& newVert) {
	for (unsigned int i = 0; i < adjacentVertices.size(); i++) {
		Vertex* temp = adjacentVertices[i];
		if ((temp)->position[0] == oldVert->position[0] &&
			(temp)->position[1] == oldVert->position[1] &&
			(temp)->position[2] == oldVert->position[2]) {
			adjacentVertices[i] = newVert;
			return;
		}
	}
}

void Vertex::collapsed(Vertex*& mid) {
	replacedBy = mid;
	for (unsigned int i = 0; i < adjacentVertices.size(); i++) {
		(adjacentVertices[i])->setNewVertForVert(this, mid);
	}
	for (unsigned int i = 0; i < adjacentFaces.size(); i++) {
		(adjacentFaces[i])->setNewVertForFace(this, mid);
	}
}

void Vertex::revertVertex() {
	isUsed = true;
	for (int i = 0; i < adjacentVertices.size(); i++) {
		for (int k = 0; k < adjacentVertices[i]->adjacentVertices.size(); k++) {
			if (adjacentVertices[i]->adjacentVertices[k]->index == replacedBy->index) {
				adjacentVertices[i]->adjacentVertices[k]->setNewVertForVert(replacedBy, this);
				break;
			}
		}
	}
	for (int j = 0; j < adjacentFaces.size(); j++) {
			adjacentFaces[j]->setNewVertForFace(replacedBy, this);
			if (adjacentFaces[j]->isDegenerate() == false) {
				adjacentFaces[j]->isUsed = true;
			}
	}
	replacedBy->isUsed = false;
}

//remove repeated vertices and faces in the adjacency lists and remove degenerate faces
void Vertex::cleanMesh() {
	std::sort(adjacentFaces.begin(), adjacentFaces.end());
	auto last = std::unique(adjacentFaces.begin(), adjacentFaces.end());
	adjacentFaces.erase(last, adjacentFaces.end());

	std::sort(adjacentVertices.begin(), adjacentVertices.end());
	auto last2 = std::unique(adjacentVertices.begin(), adjacentVertices.end());
	adjacentVertices.erase(last2, adjacentVertices.end());

	removeDegenerate();

	updateVertex();
}

void Vertex::removeDegenerate() {
	for (unsigned int i = 0; i < adjacentFaces.size(); i++) {
		Face* temp = adjacentFaces[i];
		if (temp->isUsed == false) {
			continue;
		}
		if ((temp)->isDegenerate() == true) {
			adjacentFaces[i]->isUsed = false;
		}
	}
}

void Vertex::setNormal() {
	vec3 adjacentNormalSum;
	for (unsigned int j = 0; j < adjacentFaces.size(); j++) {
		if (adjacentFaces[j]->isUsed == false) continue;
		Face* temp = adjacentFaces[j];
		adjacentNormalSum += (temp)->normal;
	}
	normal = adjacentNormalSum;
}

void Vertex::updateVertex() {
	Q = mat4x4(0);
	for (int i = 0; i < adjacentFaces.size(); i++) {
		if (adjacentFaces[i]->isUsed == false) continue;
		Q = Q + adjacentFaces[i]->K;
	}
	setNormal();
}




