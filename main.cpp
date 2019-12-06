#include "GL/glew.h"
#include "GLFW/glfw3.h"
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>
#include <tuple>
#include <complex>
//this is a test


#define M_PI 3.141592654f

unsigned int g_windowWidth = 1000;
unsigned int g_windowHeight = 800;
float modelWidth = 1.0;
float modelHeight = ((float) g_windowHeight) / ((float) g_windowWidth);

char* g_windowName = "Spin to Draw";

GLFWwindow* g_window;

// Model data
std::vector<float> g_meshVertices;
std::vector<float> g_meshNormals;
std::vector<unsigned int> g_meshIndices;
GLfloat g_modelViewMatrix[16];

// Default options
bool enablePersp = false;
bool teapotSpin = false;
bool enableDolly = false;
bool showCheckerboard = false;

// Dolly zoom options 
float fov = M_PI / 2.f;

// line up distance so that it matches the screen

float initialDistance = (1.0 / tan(fov/2)) * 0.25;
float distance = initialDistance;
// per second
float rotation_rate = 0.025f;


// these define the coordinate system
// the model is stored centered around 0, 0 and stretching from -0.5 to 0.5 along both x and y axis
std::complex<float> model_to_screen(std::complex<float> modelCoordinate) {
  return std::complex<float>(modelCoordinate.real(),
			     modelCoordinate.imag());
}

std::complex<float> mouse_to_model(double mouseX, double mouseY) {
  return std::complex<float>((mouseX - g_windowWidth/2) / g_windowWidth * modelWidth,
			     -((mouseY - g_windowHeight/2) / g_windowHeight)*modelHeight); // y axis gets flipped and center on zero
}


struct Rotating {
  std::complex<float> coefficient;
  int n;
};


struct State {
  bool mousePressed;
  std::vector<std::complex<float>> samples;
  std::vector<Rotating> rotatings;

  State() {
    mousePressed = false;
  }
};

State state;


bool compareRotating(Rotating& r1, Rotating& r2) {
  return std::abs(r1.coefficient) > std::abs(r2.coefficient);
}


const std::complex<float> complex_i(0, 1);
int granularity = 20;
int start = -(granularity / 2);

std::vector<Rotating> make_rotating(std::vector<std::complex<float>>& samples) {
  std::vector<Rotating> rotatings;
  for(int i = 0; i < granularity; i++) {
    int n = i + start;
    std::complex<float> c(0, 0);
    for(int samplei = 0; samplei < samples.size(); samplei++) {
      float t = ((float) samplei) / ((float) samples.size());
      c += samples[samplei] * std::exp(-n * M_PI * 2 * complex_i * t) / ((float) samples.size());
    }
    Rotating r;
    r.coefficient = c;
    r.n = n;
    rotatings.push_back(r);
  }

  std::sort(rotatings.begin(), rotatings.end(), compareRotating);
  
  return rotatings;
}



// Auxiliary math functions
float dotProduct(const float* a, const float* b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

void crossProduct(const float* a, const float* b, float* r)
{
	r[0] = a[1] * b[2] - a[2] * b[1];
	r[1] = a[2] * b[0] - a[0] * b[2];
	r[2] = a[0] * b[1] - a[1] * b[0];
}

float radianToDegree(float angleInRadian) {
	return angleInRadian * 180.f / M_PI;
}

void normalize(float* a)
{
	const float len = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);

	a[0] /= len;
	a[1] /= len;
	a[2] /= len;
}

void meshIndicesToVec(float* dest, int index1, int index2) {
	dest[0] = g_meshVertices[index2 * 3] - g_meshVertices[index1 * 3];
	dest[1] = g_meshVertices[index2 * 3 + 1] - g_meshVertices[index1 * 3 + 1];
	dest[2] = g_meshVertices[index2 * 3 + 2] - g_meshVertices[index1 * 3 + 2];
}

void getCurrentPosOfMouse(double &xpos, double &ypos) {
	//double xpos, ypos;
	
	glfwGetCursorPos(g_window,&xpos, &ypos);
	//std::cout << xpos << ypos << std::endl;
}

/*void computeNormals()
{
	g_meshNormals.resize(g_meshVertices.size());
	for (int i = 0; i < g_meshNormals.size(); i++) {
		g_meshNormals[i] = 0.0;
	}


	for (int v = 0; v < g_meshIndices.size() / 3; v++) {
		unsigned int index1 = g_meshIndices[v * 3];
		unsigned int index2 = g_meshIndices[v * 3 + 1];
		unsigned int index3 = g_meshIndices[v * 3 + 2];


		float vec1[3];
		float vec2[3];
		meshIndicesToVec(vec1, index1, index2);
		meshIndicesToVec(vec2, index1, index3);

		float normal[3];
		crossProduct(vec1, vec2, normal);

		g_meshNormals[index1 * 3] += normal[0];
		g_meshNormals[index1 * 3 + 1] += normal[1];
		g_meshNormals[index1 * 3 + 2] += normal[2];

		g_meshNormals[index2 * 3] += normal[0];
		g_meshNormals[index2 * 3 + 1] += normal[1];
		g_meshNormals[index2 * 3 + 2] += normal[2];

		g_meshNormals[index3 * 3] += normal[0];
		g_meshNormals[index3 * 3 + 1] += normal[1];
		g_meshNormals[index3 * 3 + 2] += normal[2];

	}

	float* vecData = g_meshNormals.data();

	for (int v = 0; v < g_meshNormals.size() / 3; v++) {
		normalize(vecData + (v * 3));

	}
}*/

double getTime()
{
	return glfwGetTime();
}

void glfwErrorCallback(int error, const char* description)
{
	std::cerr << "GLFW Error " << error << ": " << description << std::endl;
	exit(1);
}

void togglePerspective() {

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	float halfWidth = distance * tan(fov / 2);

	// Perspective Projection
	if (enablePersp)
	{

		float fovInDegree = radianToDegree(fov);
		gluPerspective(fovInDegree, (GLfloat)g_windowWidth / (GLfloat)g_windowHeight, 1.0f, 40.f);
	}
	// Othogonal Projection
	else
	{
		// Scale down the object for a better view in orthographic projection
		glScalef(0.5, 0.5, 0.5);
		//std::cout << halfWidth;

		glOrtho(-halfWidth, halfWidth, -halfWidth * g_windowHeight / g_windowWidth, halfWidth * g_windowHeight / g_windowWidth, 1.0f, 40.f);
	}
}


void mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
  if(button == GLFW_MOUSE_BUTTON_LEFT) 
    {
      // reset state on button press
      if(action == GLFW_PRESS) {
	state = State();
      } else {
	// calculate vectors on release
	state.rotatings = make_rotating(state.samples);
      }

      state.mousePressed = (action == GLFW_PRESS);
    }
}

void glfwKeyCallback(GLFWwindow* p_window, int p_key, int p_scancode, int p_action, int p_mods)
{
  
  if (p_key == GLFW_KEY_ESCAPE && p_action == GLFW_PRESS)
    {
      glfwSetWindowShouldClose(g_window, GL_TRUE);
    }
  if (p_key == GLFW_KEY_P && p_action == GLFW_PRESS) {

    // Perspective Projection
    enablePersp = true;
    togglePerspective();
    std::cout << "Perspective activated\n";

  }
  if (p_key == GLFW_KEY_O && p_action == GLFW_PRESS) {

    // Orthographic Projection
    enablePersp = false;
    togglePerspective();
    std::cout << "Orthogonal activated\n";

  }
  if (p_key == GLFW_KEY_S && p_action == GLFW_PRESS) {

    // Toggle Spinning
    if (!teapotSpin) {
      std::cout << "Teapot spinning on\n";
    }
    else {
      std::cout << "Teapot spinning off\n";
    }
    teapotSpin = !teapotSpin;
  }
  if (p_key == GLFW_KEY_D && p_action == GLFW_PRESS) {

    // Toggle dolly zoom
    if (!enableDolly)
      {
	std::cout << "Dolly zoom on\n";
      }
    else {
      std::cout << "Dolly zoom off\n";
    }
    enableDolly = !enableDolly;
  }
  if (p_key == GLFW_KEY_C && p_action == GLFW_PRESS) {

    // Show/hide Checkerboard
    if (!showCheckerboard)
      {
	std::cout << "Show checkerboard\n";
      }
    else {
      std::cout << "Hide checkerboard\n";
    }
    showCheckerboard = !showCheckerboard;
  }
}
void initWindow()
{
  // initialize GLFW
  glfwSetErrorCallback(glfwErrorCallback);
  if (!glfwInit())
    {
      std::cerr << "GLFW Error: Could not initialize GLFW library" << std::endl;
      exit(1);
    }

  g_window = glfwCreateWindow(g_windowWidth, g_windowHeight, g_windowName, NULL, NULL);
  if (!g_window)
    {
      glfwTerminate();
      std::cerr << "GLFW Error: Could not initialize window" << std::endl;
      exit(1);
    }


  // callbacks
  glfwSetKeyCallback(g_window, glfwKeyCallback);
  glfwSetMouseButtonCallback(g_window, mouse_button_callback);

  // Make the window's context current
  glfwMakeContextCurrent(g_window);

  // turn on VSYNC
  glfwSwapInterval(1);
}

void initGL()
{
  glClearColor(1.f, 1.f, 1.f, 1.0f);

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_DEPTH_TEST);
  glShadeModel(GL_SMOOTH);
}
void clearModelViewMatrix()
{
	for (int i = 0; i < 4; ++i)
	{
		for (int j = 0; j < 4; ++j)
		{
			g_modelViewMatrix[4 * i + j] = 0.0f;
		}
	}
}

void updateModelViewMatrix()
{
	clearModelViewMatrix();


	// You can use getTime() to change rotation over time
	float rotation = getTime() * 2 * M_PI * rotation_rate;

	g_modelViewMatrix[0] = 1.0f;
	g_modelViewMatrix[5] = 1.0f;
	g_modelViewMatrix[10] = 1.0f;

	if (teapotSpin) {
		// rotation matrix- keep the y axis
		g_modelViewMatrix[0] = cos(rotation);
		g_modelViewMatrix[2] = -sin(rotation);
		g_modelViewMatrix[8] = -g_modelViewMatrix[2];
		g_modelViewMatrix[10] = g_modelViewMatrix[0];
	}


	g_modelViewMatrix[14] = -distance;
	g_modelViewMatrix[15] = 1.0f;
}

void setModelViewMatrix()
{
	glMatrixMode(GL_MODELVIEW);
	updateModelViewMatrix();
	glLoadMatrixf(g_modelViewMatrix);
}

//this was a conflic and I'm not sure what the deal is so I'll comment it out for now
/*const std::complex<float> complex_i(0, 1);
int granularity = 11;
int start = -(granularity / 2);

std::vector<Rotating> make_rotating(std::vector<std::complex<float>>& samples) {
  std::vector<Rotating> rotatings;
  for(int i = 0; i < granularity; i++) {
    int n = i + start;
    std::complex<float> c(0, 0);
    for(int samplei = 0; samplei < samples.size(); samplei++) {
      float t = ((float) samplei) / ((float) samples.size());
      c += samples[samplei] * std::exp(-n * M_PI * 2 * complex_i * t) / ((float) samples.size());
    }
    Rotating r;
    r.coefficient = c;
    r.n = n;
    rotatings.push_back(r);
  }

  std::sort(rotatings.begin(), rotatings.end(), compareRotating);
  
  return rotatings;
}
std::vector <std::complex<float>> functionReturningSamples(){
	std::vector<double> vectorOfNumbers = { 2.67,1.68,-0.39,-0.26,-0.8,-0.32,-0.97,-0.15,-0.23,0.22,-0.03,0.81,0,0.9,0.11,0.32,0.26,0.4,0.21,0.55
,1.83, 3.25, 1.36, 3.23, 1.33, 3.45, -0.02, 0.16, 0.19, 0.32, 0.3, 0.4,1.88, 4.04, 2.16, 4.06, 2.3, 4.07,0.15, 0.01, 0.48, 0.03, 0.73, -0.13
,3.39, 3.7, 3.53, 3.1, 3.35, 2.48,3.6, 1.72,3.3, 1.75,2.98, 2.09,2.59, 2.32,2.25l,0.41,-1.02,2.67, 1.68 };
	std::vector<std::complex<float>> samples;
	for (int i = 0; i < vectorOfNumbers.size() - 2; i += 2) {
		samples.push_back(std::complex<float>(vectorOfNumbers[i], vectorOfNumbers[i + 1]));
	}
	samples.push_back(std::complex<float>(vectorOfNumbers[0], vectorOfNumbers[1]));
	return samples;
}

void renderLines() {
	double  xpos, ypos;
	getCurrentPosOfMouse(xpos, ypos);
	std::complex<float> newPos(xpos, ypos);
  //std::vector<std::complex<float>> samples = {std::complex<float>(-2.0, 2.0), std::complex<float>(0, 2.0),
		//				std::complex<float>(2.0, 2.0),
			//	std::complex<float>(2.0, 0.0), std::complex<float>(2.0, -2.0),
				//std::complex<float>(0.0, -2.0), std::complex<float>(-2.0, -2.0)};
	std::vector<std::complex<float>> samples = functionReturningSamples();

  std::vector<Rotating> rotatings = make_rotating(samples);*/
  
void renderLines() {
  

  glDisable(GL_LIGHTING);

  glColor3f(0, 0, 0);
  glLineWidth(3);
  glBegin(GL_LINES);
  
  float time = getTime() * 2 * M_PI * rotation_rate;
  std::complex<float> currentpos(0,0);

  for(Rotating r : state.rotatings) {
    int n = r.n;
    std::complex<float> oldpos(currentpos);
    currentpos += r.coefficient * std::exp(n * 2 * M_PI * complex_i * time);

    std::complex<float> oldposScreen = model_to_screen(oldpos);
    std::complex<float> currentposScreen = model_to_screen(currentpos);
    glVertex3f(oldposScreen.real(), oldposScreen.imag(), 0);
    glVertex3f(currentposScreen.real(), currentposScreen.imag(), 0);
  }

  glEnd();

  glEnable(GL_LIGHTING);
}


void renderDrawing() {
  glDisable(GL_LIGHTING);

  glColor3f(0, 0, 0);
  glLineWidth(3);
  glBegin(GL_LINES);
  

  for(int i = 0; i < state.samples.size(); i++) {
    int lastIndex = i - 1;
    if(lastIndex < 0) {
      lastIndex = state.samples.size() -1;
    }
    
    std::complex<float> last = model_to_screen(state.samples[lastIndex]);
    std::complex<float> current = model_to_screen(state.samples[i]);
    glVertex3f(last.real(), last.imag(), 0);
    glVertex3f(current.real(), current.imag(), 0);
  }
  
  glEnd();

  glEnable(GL_LIGHTING);
}

void drawCheckerBoard() {
	int checkerCount = g_windowWidth;
	int checkerSize = (g_windowWidth) / checkerCount;

	glBegin(GL_QUADS);
	for (int i = 0; i < checkerCount; ++i) {
		for (int j = 0; j < checkerCount; ++j) {
			if ((i + j) % 2 == 0)
				glColor3f(0.0, 0.1, 1.0);
			else
				glColor3f(1.0, 0.0, 1.0);

			float x = i - checkerCount / 2; // to be between -1 and 1
			float z = j - checkerCount / 2;
			x *= checkerSize;
			z *= checkerSize;
			float y = -1.0f;
			glVertex3f(x, y, z);
			glVertex3f(x, y, z - checkerSize);
			glVertex3f(x + checkerSize, y, z - checkerSize);
			glVertex3f(x + checkerSize, y, z);
		}
	}
	glEnd();
}

void renderCheckerBoard() {

	/*
	/* If you want to keep checkerboard still while rotating the
	/* the teapot, you need to change the transformation for the
	/* checkerboard plane
	*/
	glMatrixMode(GL_MODELVIEW);
	clearModelViewMatrix();

	g_modelViewMatrix[0] = 1;
	g_modelViewMatrix[2] = 0;
	g_modelViewMatrix[5] = 1;
	g_modelViewMatrix[8] = 0;
	g_modelViewMatrix[10] = 1;
	g_modelViewMatrix[14] = -distance;
	g_modelViewMatrix[15] = 1.0f;

	glLoadMatrixf(g_modelViewMatrix);

	// Disable shading for the checkerboard
	glDisable(GL_LIGHTING);
	drawCheckerBoard();
	glEnable(GL_LIGHTING);
}

void render()
{
	togglePerspective();
	setModelViewMatrix();

	renderDrawing();
	if (!state.mousePressed) {
	  renderLines();
	}
	

	if (showCheckerboard)
		renderCheckerBoard();
}

void onTick(GLFWwindow* window) {
  if (state.mousePressed) {
    double xpos, ypos;
    glfwGetCursorPos(window, &xpos, &ypos);
    
    state.samples.push_back(mouse_to_model(xpos, ypos));
  }
}


void renderLoop()
{
	while (!glfwWindowShouldClose(g_window))
	{
		// clear buffers
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		render();

		// Swap front and back buffers
		glfwSwapBuffers(g_window);

		// Poll for and process events
		glfwPollEvents();
		
		onTick(g_window);
	}
}

int main()
{
	initWindow();
	initGL();
	printHotKeys();
	renderLoop();
}
