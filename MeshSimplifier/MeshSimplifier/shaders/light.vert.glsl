# version 330 core

layout (location = 0) in vec3 position;
layout (location = 1) in vec3 normal;
layout (location = 2) in vec3 color;

uniform mat4 projection;
uniform mat4 modelview;

flat out vec3 Normal; 
out vec4 Position;
flat out vec4 Color;

void main() {
    gl_Position = projection * modelview * vec4(position, 1.0f);
    Normal = mat3(transpose(inverse(modelview))) * normal; 
    Position = modelview * vec4(position, 1.0f);
	Color = vec4(color, 1.0f);
}

