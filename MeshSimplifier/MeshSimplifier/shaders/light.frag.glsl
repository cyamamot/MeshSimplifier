# version 330 core

flat in vec3 Normal; 
in vec4 Position;
flat in vec4 Color;

// The output of the fragment shader will determine the fragment's color
out vec4 fragColor;

uniform int islight ; // are we lighting. 

uniform vec4 color; // RGB color normally used in glColor*

uniform vec4 light0posn ; 
uniform vec4 light0color ; 
uniform vec4 light1posn ; 
uniform vec4 light1color ; 
uniform vec4 light2posn ; 
uniform vec4 light2color ; 


uniform vec4 ambient ; 
uniform vec4 diffuse ; 
uniform vec4 specular ; 
uniform float shininess ; 

vec4 ComputeLight (vec3 direction, vec4 lightcolor, vec3 normal, vec3 halfvec, vec4 mydiffuse, vec4 myspecular, float myshininess) {

        float nDotL = dot(normal, direction)  ;         
        vec4 lambert = mydiffuse * lightcolor * max (nDotL, 0.0) ;  

        float nDotH = dot(normal, halfvec) ; 
        vec4 phong = myspecular * lightcolor * pow (max(nDotH, 0.0), myshininess) ; 

        vec4 retval = lambert + phong ; 
        return retval ;
}       

void main (void) 
{       
    if (islight == 0){
		//fragColor = vec4(color, 1.0f); 
		fragColor = Color;
	}
    else { 
        // The eye is always at (0,0,0) looking down -z axis 

        const vec3 eyepos = vec3(0,0,0) ; 
        vec3 mypos = Position.xyz / Position.w ; // Dehomogenize current location 
        vec3 eyedirn = normalize(eyepos - mypos) ; 

        // Compute normal, needed for shading. 
        vec3 normal = normalize(Normal) ; 

        // Light 0, point
        vec3 position0 = light0posn.xyz / light0posn.w ; 
        vec3 direction0 = normalize (position0 - mypos) ; // no attenuation 
        vec3 half0 = normalize (direction0 + eyedirn) ;  
        vec4 col0 = ComputeLight(direction0, light0color, normal, half0, diffuse, specular, shininess) ;

        // Light 1, point 
        vec3 position1 = light1posn.xyz / light1posn.w ; 
        vec3 direction1 = normalize (position1 - mypos) ; // no attenuation 
        vec3 half1 = normalize (direction1 + eyedirn) ;  
        vec4 col1 = ComputeLight(direction1, light1color, normal, half1, diffuse, specular, shininess) ;

		
        // Light 2, point 
        vec3 position2 = light2posn.xyz / light2posn.w ; 
        vec3 direction2 = normalize (position2 - mypos) ; // no attenuation 
        vec3 half2 = normalize (direction2 + eyedirn) ;  
        vec4 col2 = ComputeLight(direction2, light2color, normal, half2, diffuse, specular, shininess) ;
        
        fragColor = ambient + col0 + col1 + col2; 
	}
}
