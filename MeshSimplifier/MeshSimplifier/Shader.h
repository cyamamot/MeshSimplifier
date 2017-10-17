/*
Shader class provided by UCSD CSE167 courtesy Prof. Ravi Ramamoorthi
*/
#pragma once
#include <iostream>
#include <string>
#include "Dependencies/glm/glm.hpp"
#include "Dependencies/glew/glew.h"
#include "Dependencies/glut/glut.h"

#ifndef __INCLUDESHADERS 
#define __INCLUDESHADERS 

std::string textFileRead(const char * filename);
void programerrors(const GLint program);
void shadererrors(const GLint shader);
GLuint initshaders(GLenum type, const char * filename);
GLuint initprogram(GLuint vertexshader, GLuint fragmentshader);

#endif 
