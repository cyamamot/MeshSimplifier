/*
Shader class used to open and write to shaders
*/
#pragma once
#include <iostream>
#include <string>
#include "Dependencies/glm/glm.hpp"
#include "Dependencies/glew/glew.h"
#include "Dependencies/glut/glut.h"

#ifndef __INCLUDESHADERS 
#define __INCLUDESHADERS 

std::string ReadFile(const char * filename);
void LogProgramError(const GLint program);
void LogShaderError(const GLint shader);
GLuint InitializeShaders(GLenum type, const char * filename);
GLuint InitializeProgram(GLuint vertexshader, GLuint fragmentshader);

#endif 
