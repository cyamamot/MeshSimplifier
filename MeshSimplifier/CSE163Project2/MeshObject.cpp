#include "MeshObject.h"
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>  
#include <queue>


MeshObject::MeshObject(const char* filename) {
	//Constructor parses OFF format object file and creates each face and vertex
	std::vector<float> v;
	std::ifstream infile(filename);
	setVertices(v, infile);

	//centers object on screen by adjusting each vertex
	float avgY = (Ymin + Ymax) / 2.0f - 0.02f;
	float avgZ = (Zmin + Zmax) / 2.0f;
	for (unsigned int i = 0; i < vertexPositionList.size(); ++i) {
		vec3 shiftedVertex = (vertexPositionList[i] - vec3(0.0f, avgY, avgZ)) * vec3(1.58f, 1.58f, 1.58f);
		vertexPositionList[i] = shiftedVertex;
		vertList[i]->position = shiftedVertex;
	}
	v.clear();
	//create each face with 3 vertices from the list and sets adjacencies of every vertex
	setFaces(v, infile);
	//removes repeated elements from the adjacency lists of each vertex
	for (int i = 0; i < numVertices; i++) {
		vertList[i]->cleanMesh();
	}
	setVertexNormals();


	std::cout << "here" << std::endl;
}

void MeshObject::setVertices(std::vector<float>& v, std::ifstream& infile) {
	//parse first two lines (in .off files they are data on the file, not vertex data)
	if (!infile) return;
	std::string line;
	std::getline(infile, line);
	std::getline(infile, line);
	std::istringstream iss(line);
	float n;
	while (iss >> n) {
		v.push_back(n);
	}
	//creates pointer to every vertex
	numVertices = (int)v[0];
	numFaces = (int)v[1];
	v.clear();
	for (int i = 0; i < numVertices; i++) {
		std::getline(infile, line);
		std::istringstream iss(line);
		while (iss >> n) {
			v.push_back(n);
		}
		if (v[0] < Xmin) Xmin = v[0]; if (v[0] > Xmax) Xmax = v[0];
		if (v[1] < Ymin) Ymin = v[1]; if (v[1] > Ymax) Ymax = v[1];
		if (v[2] < Zmin) Zmin = v[2]; if (v[2] > Zmax) Zmax = v[2];
		Vertex *vert = new Vertex(v[0], v[1], v[2]);
		vert->index = i;
		vertList.push_back(vert);
		vertexPositionList.push_back(vec3(v[0], v[1], v[2]));
		faceColorList.push_back(vec3(1, 1, 0));
		v.clear();
	}
}

void MeshObject::setFaces(std::vector<float>& v, std::ifstream& infile) {
	std::string line;
	float n;
	for (int i = 0; i < numFaces; i++) {
		std::getline(infile, line);
		std::istringstream iss(line);
		while (iss >> n) {
			v.push_back(n);
		}
		vertexIndexList.push_back(v[1]);
		vertexIndexList.push_back(v[2]);
		vertexIndexList.push_back(v[3]);
		//if ()
		Face* face = new Face(vertList[v[1]], vertList[v[2]], vertList[v[3]]);
		face->index = i;
		faceList.push_back(face);
		vertList[v[1]]->adjacentFaces.push_back(faceList[i]);
		vertList[v[2]]->adjacentFaces.push_back(faceList[i]);
		vertList[v[3]]->adjacentFaces.push_back(faceList[i]);
		vertList[v[1]]->adjacentVertices.push_back(vertList[v[2]]);
		vertList[v[1]]->adjacentVertices.push_back(vertList[v[3]]);
		vertList[v[2]]->adjacentVertices.push_back(vertList[v[1]]);
		vertList[v[2]]->adjacentVertices.push_back(vertList[v[3]]);
		vertList[v[3]]->adjacentVertices.push_back(vertList[v[1]]);
		vertList[v[3]]->adjacentVertices.push_back(vertList[v[2]]);
		face->Color = vec3(1, 1, 0);
		v.clear();
	}
}

int MeshObject::getNumFaces() {
	return vertList.size();
}

bool areEqual(Face*& a, Face*& b) {
	if (!a->isDegenerate() && !b->isDegenerate()) {
		if (a->vertices[0]->index == b->vertices[0]->index || a->vertices[0]->index == b->vertices[1]->index ||
			a->vertices[0]->index == b->vertices[2]->index) {
			if (a->vertices[1]->index == b->vertices[0]->index || a->vertices[1]->index == b->vertices[1]->index ||
				a->vertices[1]->index == b->vertices[2]->index) {
				if (a->vertices[2]->index == b->vertices[0]->index || a->vertices[2]->index == b->vertices[1]->index ||
					a->vertices[2]->index == b->vertices[2]->index) {
					return true;
				}
			}
		}
	}
	return false;
}

//used for scaling object in scene
float MeshObject::norm() {
	float maxX = std::max(Xmax, std::abs(Xmin));
	float maxY = std::max(Ymax, std::abs(Ymin));
	float maxZ = std::max(Zmax, std::abs(Zmin));
	return std::max(maxX, std::max(maxY, maxZ));
}

void MeshObject::initialize(GLuint &VAO, GLuint &VBO, GLuint &EBO, GLuint &NBO, GLuint &CBO) {

	glBindVertexArray(VAO);

	// Bind vertices to layout location 0
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vec3) * vertexPositionList.size(), &vertexPositionList[0], GL_DYNAMIC_DRAW);
	glEnableVertexAttribArray(0); // This allows usage of layout location 0 in the vertex shader
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), 0);

	// Bind normals to layout location 1
	glBindBuffer(GL_ARRAY_BUFFER, NBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vec3) * normalsList.size(), &normalsList[0], GL_DYNAMIC_DRAW);
	glEnableVertexAttribArray(1); // This allows usage of layout location 1 in the vertex shader
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), 0);

	// Bind colors to layout location 2
	glBindBuffer(GL_ARRAY_BUFFER, CBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vec3) * faceColorList.size(), &faceColorList[0], GL_DYNAMIC_DRAW);
	glEnableVertexAttribArray(2); // This allows usage of layout location 2 in the vertex shader
	glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), 0);

	// The indices array tells OpenGL what order to iterate through the buffers in when the shaders execute
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(int) * vertexIndexList.size(), &vertexIndexList[0], GL_DYNAMIC_DRAW);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);
}

void MeshObject::drawObject(GLuint &VAO, GLuint &VBO, GLuint &EBO, GLuint &NBO, GLuint &CBO) {
	glUniformMatrix4fv(modelviewPos, 1, GL_FALSE, &(modelview)[0][0]);
	initialize(VAO, VBO, EBO, NBO, CBO);
	glBindVertexArray(VAO);
	glDrawElements(GL_TRIANGLES, vertexIndexList.size(), GL_UNSIGNED_INT, 0);
	glBindVertexArray(0); // Unbind the VAO when done
}

//Helper method to set the normal of each vertex based on the normals of the adjacent faces
void MeshObject::setVertexNormals() {
	for (int i = 0; i < vertList.size(); i++) {
		/*if (vertList[i]->isUsed == false) {
			normalsList.push_back(vec3(0, 0, 0));
			continue;
		}*/
		vec3 adjacentNormalSum;

		for (unsigned int j = 0; j < vertList[i]->adjacentFaces.size(); j++) {
			if (vertList[i]->adjacentFaces[j]->isUsed == false) continue;
			Face* temp = vertList[i]->adjacentFaces[j];
			adjacentNormalSum += (temp)->normal;
		}
		vertList[i]->normal = adjacentNormalSum * (float)(1.0f / vertList[i]->adjacentFaces.size());
		normalsList.push_back(-vertList[i]->normal);
	}
}


//collapses an edge between two specified vertices
void MeshObject::edgeCollapse(int vert1, int vert2) {
	if (vert1 >= vertList.size() || vert1 < 0 || vert2 >= vertList.size() || vert2 < 0) {
		std::cout << "vertices with these indices are not present" << std::endl;
		return;
	}
	Vertex *vertex1 = vertList[vert1];
	Vertex *vertex2 = vertList[vert2];
	if (vertex1->isUsed == false || vertex2->isUsed == false) {
		std::cout << "vertices not in use" << std::endl;
	}
	bool verticesConnected = false;
	for (int i = 0; i < vertex1->adjacentVertices.size(); i++) {
		if (vertex1->adjacentVertices[i]->index == vert2) verticesConnected = true;
	}
	if (verticesConnected == false) {
		std::cout << "vertices are not connected" << std::endl;
		return;
	}
	double x = vertex1->position[0] + vertex2->position[0];
	double y = vertex1->position[1] + vertex2->position[1];
	double z = vertex1->position[2] + vertex2->position[2];
	//create new vertex between two collapsed vertices
	Vertex *midpoint = new Vertex(x / 2.0, y / 2.0, z / 2.0);

	midpoint->absorb(vertex1, vertex2);
	midpoint->index = vertList.size();
	vertList.push_back(midpoint);
	vertex1->collapsed(midpoint);
	vertex1->replacedBy = midpoint;
	vertList[vert1]->isUsed = false;
	vertex2->collapsed(midpoint);
	vertex2->replacedBy = midpoint;
	vertList[vert2]->isUsed = false;
	vertexPositionList.push_back(midpoint->position);
	faceColorList.push_back(vec3(1, 1, 0));
	midpoint->cleanMesh();
	for (int i = 0; i < midpoint->adjacentFaces.size(); i++) {
		for (int j = i + 1; j < midpoint->adjacentFaces.size(); j++) {
			if (midpoint->adjacentFaces[i]->isUsed && midpoint->adjacentFaces[j]->isUsed &&
						areEqual(midpoint->adjacentFaces[i], midpoint->adjacentFaces[j])) {
				midpoint->adjacentFaces[i]->setNotUsed();
				midpoint->adjacentFaces[j]->setNotUsed();
			}
		}
	}
	normalsList.push_back(-midpoint->normal);
	for (int i = 0; i < midpoint->adjacentVertices.size(); i++) {
		if (midpoint->adjacentVertices[i]->isUsed) {
			midpoint->adjacentVertices[i]->cleanMesh();
			normalsList[midpoint->adjacentVertices[i]->index] = -midpoint->adjacentVertices[i]->normal;
		}
	}
	for (int i = 0; i < midpoint->adjacentFaces.size(); i++) {
		if (midpoint->adjacentFaces[i]->isUsed == false) {
			vertexIndexList[3 * midpoint->adjacentFaces[i]->index] = -1;
			vertexIndexList[3 * midpoint->adjacentFaces[i]->index + 1] = -1;
			vertexIndexList[3 * midpoint->adjacentFaces[i]->index + 2] = -1;
			continue;
		}
		vertexIndexList[3 * midpoint->adjacentFaces[i]->index] = midpoint->adjacentFaces[i]->vertices[0]->index;
		vertexIndexList[3 * midpoint->adjacentFaces[i]->index + 1] = midpoint->adjacentFaces[i]->vertices[1]->index;
		vertexIndexList[3 * midpoint->adjacentFaces[i]->index + 2] = midpoint->adjacentFaces[i]->vertices[2]->index;
	}
}

//undoes a collapsed edge
void MeshObject::undoCollapse(int vert1, int vert2) {
	if (vert1 >= vertList.size() || vert1 < 0 || vert2 >= vertList.size() || vert2 < 0) {
		std::cout << "vertices with these indices are not present" << std::endl;
		return;
	}
	Vertex *vertex1 = vertList[vert1];
	Vertex *vertex2 = vertList[vert2];
	if (vertex1->isUsed || vertex2->isUsed) return;
	vertex1->revertVertex();
	vertex2->revertVertex();

	vertList.pop_back();
	vertexPositionList.pop_back();
	faceColorList.pop_back();
	normalsList.pop_back();

	for (int i = 0; i < vertex1->adjacentVertices.size(); i++) {
		if (vertex1->adjacentVertices[i]->isUsed) {
			vertex1->adjacentVertices[i]->cleanMesh();
			normalsList[vertex1->adjacentVertices[i]->index] = -vertex1->adjacentVertices[i]->normal;
		}
	}
	for (int i = 0; i < vertex2->adjacentVertices.size(); i++) {
		if (vertex2->adjacentVertices[i]->isUsed) {
			vertex2->adjacentVertices[i]->cleanMesh();
			normalsList[vertex2->adjacentVertices[i]->index] = -vertex2->adjacentVertices[i]->normal;
		}
	}
	for (int i = 0; i < vertex1->adjacentFaces.size(); i++) {
		if (vertex1->adjacentFaces[i] == false) {
			vertexIndexList[3 * vertex1->adjacentFaces[i]->index] = -1;
			vertexIndexList[3 * vertex1->adjacentFaces[i]->index + 1] = -1;
			vertexIndexList[3 * vertex1->adjacentFaces[i]->index + 2] = -1;
			continue;
		}
		vertexIndexList[3 * vertex1->adjacentFaces[i]->index] = vertex1->adjacentFaces[i]->vertices[0]->index;
		vertexIndexList[3 * vertex1->adjacentFaces[i]->index + 1] = vertex1->adjacentFaces[i]->vertices[1]->index;
		vertexIndexList[3 * vertex1->adjacentFaces[i]->index + 2] = vertex1->adjacentFaces[i]->vertices[2]->index;
	}
	for (int i = 0; i < vertex2->adjacentFaces.size(); i++) {
		if (vertex2->adjacentFaces[i] == false) {
			vertexIndexList[3 * vertex2->adjacentFaces[i]->index] = -1;
			vertexIndexList[3 * vertex2->adjacentFaces[i]->index + 1] = -1;
			vertexIndexList[3 * vertex2->adjacentFaces[i]->index + 2] = -1;
			continue;
		}
		vertexIndexList[3 * vertex2->adjacentFaces[i]->index] = vertex2->adjacentFaces[i]->vertices[0]->index;
		vertexIndexList[3 * vertex2->adjacentFaces[i]->index + 1] = vertex2->adjacentFaces[i]->vertices[1]->index;
		vertexIndexList[3 * vertex2->adjacentFaces[i]->index + 2] = vertex2->adjacentFaces[i]->vertices[2]->index;
	}
}

void MeshObject::updateVBO(GLuint &VAO, GLuint &VBO, GLuint &EBO, GLuint &NBO, GLuint &CBO) {
	glBindVertexArray(VAO);

	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	vec3 *mapped = reinterpret_cast<glm::vec3*>(glMapBufferRange(GL_ARRAY_BUFFER, 0, sizeof(glm::vec3)*vertexPositionList.size(),
		GL_MAP_WRITE_BIT | GL_MAP_INVALIDATE_BUFFER_BIT));
	std::copy(vertexPositionList.begin(), vertexPositionList.end(), mapped);
	glUnmapBuffer(GL_ARRAY_BUFFER);

	// Bind normals to layout location 1
	glBindBuffer(GL_ARRAY_BUFFER, NBO);
	mapped = reinterpret_cast<glm::vec3*>(glMapBufferRange(GL_ARRAY_BUFFER, 0, sizeof(glm::vec3)*normalsList.size(),
		GL_MAP_WRITE_BIT | GL_MAP_INVALIDATE_BUFFER_BIT));
	std::copy(normalsList.begin(), normalsList.end(), mapped);
	glUnmapBuffer(GL_ARRAY_BUFFER);

	// Bind colors to layout location 2
	glBindBuffer(GL_ARRAY_BUFFER, CBO);
	mapped = reinterpret_cast<glm::vec3*>(glMapBufferRange(GL_ARRAY_BUFFER, 0, sizeof(glm::vec3)*faceColorList.size(),
		GL_MAP_WRITE_BIT | GL_MAP_INVALIDATE_BUFFER_BIT));
	std::copy(faceColorList.begin(), faceColorList.end(), mapped);
	glUnmapBuffer(GL_ARRAY_BUFFER);

	// The indices array tells OpenGL what order to iterate through the buffers in when the shaders execute
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	int *mapped2 = reinterpret_cast<int*>(glMapBufferRange(GL_ELEMENT_ARRAY_BUFFER, 0, sizeof(int)*vertexIndexList.size(),
		GL_MAP_WRITE_BIT | GL_MAP_INVALIDATE_BUFFER_BIT));
	std::copy(vertexIndexList.begin(), vertexIndexList.end(), mapped2);
	glUnmapBuffer(GL_ELEMENT_ARRAY_BUFFER);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);

}

struct MeshObject::quadricContainer {
	Vertex* first;
	Vertex* second;
	float error;
	quadricContainer(Vertex*& a, Vertex*& b, float c) : first(a), second(b), error(c) {}
};

void MeshObject::printMatrix(mat4x4& m) {
	std::cout << m[0][0] << "   " << m[1][0] << "   " << m[2][0] << "   " << m[3][0] << std::endl;
	std::cout << m[0][1] << "   " << m[1][1] << "   " << m[2][1] << "   " << m[3][1] << std::endl;
	std::cout << m[0][2] << "   " << m[1][2] << "   " << m[2][2] << "   " << m[3][2] << std::endl;
	std::cout << m[0][3] << "   " << m[1][3] << "   " << m[2][3] << "   " << m[3][3] << std::endl;
	std::cout << std::endl;
}

struct MeshObject::Comp {
	bool operator()(quadricContainer*& a, quadricContainer*& b) {
		Vertex* E1V1 = a->first;
		Vertex* E1V2 = a->second;
		Vertex* E2V1 = b->first;
		Vertex* E2V2 = b->second;

		Vertex E1Mid = Vertex((E1V1->position[0] + E1V2->position[0]) / 2.0, (E1V1->position[1] + E1V2->position[1]) / 2.0, (E1V1->position[2] + E1V2->position[2]) / 2.0);
		Vertex E2Mid = Vertex((E2V1->position[0] + E2V2->position[0]) / 2.0, (E2V1->position[1] + E2V2->position[1]) / 2.0, (E2V1->position[2] + E2V2->position[2]) / 2.0);

		mat4x4 E1Q = E1V1->Q + E1V2->Q;
		mat4x4 E2Q = E2V1->Q + E2V2->Q;
		float E1Error = glm::dot((vec4(E1Mid.position, 1) * E1Q), vec4(E1Mid.position, 1));
		float E2Error = glm::dot((vec4(E2Mid.position, 1) * E2Q), vec4(E2Mid.position, 1));

		a->error = E1Error;
		b->error = E2Error;
		
		return E1Error > E2Error;
	}
};

//uses quadric simplification from Garland's '97 paper on mesh simplification 
//applies simplification method simplificationIterations(25) number of times
void MeshObject::quadricSimplification() {
	for (int a = 0; a < simplificationIterations; a++) {
		std::vector<quadricContainer*> heap;
		for (int i = 0; i < vertList.size(); i++) {
			if (vertList[i]->isUsed == false) continue;
			for (int j = 0; j < vertList[i]->adjacentVertices.size(); j++) {
				float err = 0;
				if (vertList[i]->adjacentVertices[j]->isUsed == true) heap.push_back(new quadricContainer(vertList[i], vertList[i]->adjacentVertices[j], err));
			}
		}
		if (heap.empty()) return;
		std::make_heap(heap.begin(), heap.end(), Comp());
		edgeCollapse(heap.front()->first->index, heap.front()->second->index);
		unsimplificationOrder.push_front(heap.front()->second->index);
		unsimplificationOrder.push_front(heap.front()->first->index);
	}
}

void MeshObject::quadricUnsimplification() {
	if (unsimplificationOrder.empty()) return;
	for (int a = 0; a < simplificationIterations; a++) {
		int v1 = unsimplificationOrder.front();
		unsimplificationOrder.pop_front();
		int v2 = unsimplificationOrder.front();
		unsimplificationOrder.pop_front();
		undoCollapse(v1, v2);
	}
}

