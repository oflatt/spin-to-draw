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


#define M_PI 3.141592654f

unsigned int g_windowWidth = 1000;
unsigned int g_windowHeight = 800;
float modelWidth = 1.0;
float modelHeight = ((float) g_windowHeight) / ((float) g_windowWidth);

char* g_windowName = "HW3-3D-Basics";

GLFWwindow* g_window;
std::string lastThingStopped;
// Model data
std::vector<float> g_meshVertices;
std::vector<float> g_meshNormals;
std::vector<unsigned int> g_meshIndices;
GLfloat g_modelViewMatrix[16];

// Default options
bool enablePersp = false;
bool teapotSpinLeft = false;
bool teapotSpinRight = false;
bool enableDolly = false;
bool showCheckerboard = false;

// Dolly zoom options 
float fov = M_PI / 2.f;

// line up distance so that it matches the screen

float initialDistance = (1.0 / tan(fov/2)) * 0.25;
float distance = initialDistance;

// per second
float rotation_rate = 0.025f;
void renderDrawing();

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
  double time = 0;
  float rotation;
  GLfloat g_modelViewMatrixState[16];
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
		  state.rotatings = make_rotating(state.samples);
	// calculate vectors on release
      }

      state.mousePressed = (action == GLFW_PRESS);
    }
}
void getCurrentPosOfMouse(double &xpos, double &ypos) {
	//double xpos, ypos;

	glfwGetCursorPos(g_window, &xpos, &ypos);
	//std::cout << xpos << ypos << std::endl;
}
void glfwKeyCallback(GLFWwindow* p_window, int p_key, int p_scancode, int p_action, int p_mods)
{
  
  if (p_key == GLFW_KEY_ESCAPE && p_action == GLFW_PRESS)
    {
      glfwSetWindowShouldClose(g_window, GL_TRUE);
    }
  if (p_key == GLFW_KEY_RIGHT && p_action == GLFW_PRESS) {

	  // Toggle Spinning
	  if (!teapotSpinRight) {
		  std::cout << "Teapot spinning on\n";
	  }
	  else {
		  lastThingStopped = "right";
		  std::cout << "Teapot spinning off\n";
	  }
	  teapotSpinRight = !teapotSpinRight;
	  teapotSpinLeft = false;
  }
  if (p_key == GLFW_KEY_O && p_action == GLFW_PRESS) {

    // Orthographic Projection
    enablePersp = false;
    togglePerspective();
    std::cout << "Orthogonal activated\n";

  }
  if (p_key == GLFW_KEY_LEFT && p_action == GLFW_PRESS) {

    // Toggle Spinning
    if (!teapotSpinLeft) {
      std::cout << "Teapot spinning on\n";
    }
    else {
	  lastThingStopped = "left";
      std::cout << "Teapot spinning off\n";
    }
    teapotSpinLeft = !teapotSpinLeft;
	teapotSpinRight = false;
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

void printHotKeys() {
	std::cout << "\nHot Keys..\n"
		<< "Orthogonal Projection:  O\n"
		<< "Perspective Projection: P\n"
		<< "Toggle Spinning:        S\n"
		<< "Toggle Dolly Zoom:      D\n"
		<< "Show/hide Checkerboard: C\n"
		<< "Exit:                   Esc\n\n";
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

	//state.time = abs(state.time - getTime());
	float rotation = state.rotation;
	if (lastThingStopped == "right") {

		g_modelViewMatrix[0] = cos(1 * rotation);
		g_modelViewMatrix[2] = sin(1 * rotation);
		g_modelViewMatrix[5] = 1;
		g_modelViewMatrix[8] = sin(1 * rotation);
		g_modelViewMatrix[10] = -cos(1 * rotation);
	}
	else {
		g_modelViewMatrix[0] = cos(1 * rotation);
		g_modelViewMatrix[2] = -sin(1 * rotation);
		g_modelViewMatrix[5] = 1;
		g_modelViewMatrix[8] = sin(1 * rotation);
		g_modelViewMatrix[10] = cos(1 * rotation);
	}

	if (teapotSpinRight) {
		state.rotation = state.time * 2 * M_PI * rotation_rate;
		float rotation = state.rotation;
		g_modelViewMatrix[0] = cos(1 * rotation);
		g_modelViewMatrix[2] = sin(1 * rotation);
		g_modelViewMatrix[5] = 1;
		g_modelViewMatrix[8] = sin(1 * rotation);
		g_modelViewMatrix[10] = -cos(1 * rotation);
		state.time = state.time + rotation_rate + 0.15;
		//state.time = getTime()-state.time;
	}
	else if (teapotSpinLeft) {
		state.rotation = (state.time * 2 * M_PI * rotation_rate);
		float rotation = state.rotation;
		g_modelViewMatrix[0] = cos(1 * rotation);
		g_modelViewMatrix[2] = -sin(1 * rotation);
		g_modelViewMatrix[5] = 1;
		g_modelViewMatrix[8] = sin(1 * rotation);
		g_modelViewMatrix[10] = cos(1 * rotation);
		state.time = state.time + rotation_rate + 0.15;

	}

	g_modelViewMatrix[14] = -distance;
	g_modelViewMatrix[15] = 1.0f;
	//state.time = abs(state.time - getTime());

}
void setModelViewMatrix()
{
	glMatrixMode(GL_MODELVIEW);
	updateModelViewMatrix();
	glLoadMatrixf(g_modelViewMatrix);
}


void renderLines() {
  

  glDisable(GL_LIGHTING);

  glColor3f(0, 0, 0);
  glLineWidth(3);
  glBegin(GL_LINES);
  
  float time = getTime() * 2 * M_PI * rotation_rate;
  std::complex<float> currentpos(0, 0);

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
