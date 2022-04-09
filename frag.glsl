#version 460 core

// All of the following variables could be defined in the OpenGL
// program and passed to this shader as uniform variables. This
// would be necessary if their values could change during runtim.
// However, we will not change them and therefore we define them 
// here for simplicity.

//vec3 I = vec3(1, 1, 1);          // point light intensity
vec3 I = vec3(0.6, 0.6, 0.6);          // point light intensity
vec3 Iamb = vec3(0.8, 0.8, 0.8); // ambient light intensity
vec3 kd = vec3(1, 0.2, 0.2);     // diffuse reflectance coefficient
vec3 ka = vec3(0.3, 0.3, 0.3);   // ambient reflectance coefficient
vec3 ks = vec3(0.8, 0.8, 0.8);   // specular reflectance coefficient
//vec3 lightPos = vec3(7.5, 7, 2);   // light position in world coordinates
vec3 lightPos = vec3(5, 5, 5);   // light position in world coordinates

uniform vec3 eyePos;
uniform sampler2D flagTexture;

in vec4 fragWorldPos;
in vec3 fragWorldNor;
in vec2 fragTex;

out vec4 fragColor;

void main(void)
{
	// Compute lighting. We assume lightPos and eyePos are in world
	// coordinates. fragWorldPos and fragWorldNor are the interpolated
	// coordinates by the rasterizer.

	vec3 L = normalize(lightPos - vec3(fragWorldPos));
	vec3 V = normalize(eyePos - vec3(fragWorldPos));
	vec3 H = normalize(L + V);
	vec3 N = normalize(fragWorldNor);

	float NdotL = dot(N, L); // for diffuse component
	float NdotH = dot(N, H); // for specular component

    //kd = texture(flagTexture, vec2(260.f/280.f, 60.f/178.f)).xyz;
    kd = texture(flagTexture, fragTex).xyz;

	vec3 diffuseColor = I * kd * max(0, NdotL);
	vec3 specularColor = I * ks * pow(max(0, NdotH), 100);
	vec3 ambientColor = Iamb * ka;


	fragColor = vec4(diffuseColor + specularColor + ambientColor, 1.0f);

    //fragColor = vec4(0.f,0.f, 5*kd.y, 1.f);
    //if(fragTex.s > 0.4f)
    //{
    //    fragColor = vec4(0.6f, 0.f, 0.6f, 1.f);
    //}
    //if(fragTex.t > 0.55f)
    //{
    //    fragColor = vec4(0.8f, 0.6f, 0.f, 1.f);
    //}
        //fragColor = vec4(0.6f, 0.f, 0.6f, 1.f);
    
	//fragColor = vec4(texturePixel.xyz * fragColor.xyz, 1.0f);
    //fragColor = vec4(0.8f, 0.2f, 0.2f, 1);
}
