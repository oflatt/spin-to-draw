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

const std::complex<float> complex_i(0, 1);
int granularity = 20;
int start = -(granularity / 2);



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


std::complex<float> circleFunction(float time) {
  return std::exp(M_PI * 2 * complex_i * time);
}


// lines specified by two points each
// returns a float that is the magnitude to multiply point2-point1 by
// a + bt = cj the x part of the vectors, a is the constant
// ap + bp(t) = cp(j) solved
std::pair<bool, float> interceptLines(std::complex<float> point1, std::complex<float> point2, std::complex<float> point3, std::complex<float> point4) {
  float a = point1.real() - point3.real();
  float b = point2.real()-point1.real();
  float c = point4.real();
  float ap = point1.imag() - point3.imag();
  float bp = point2.imag()-point1.imag();
  float cp = point4.imag();
  if (c == 0.0) {
	  if (b != 0.0) {
		  return std::make_pair(false, 0.0);
	  }
	  else {
		  return std::make_pair(true, -a / b);
	  }
  }
  if(bp-b*cp == 0.0) {
    return std::make_pair(false, 0.0);
  }
  float firstRes = c / (c*bp - b*cp);
  float secondRes = (cp*a)/c - ap;
  return std::make_pair(true, firstRes * secondRes);
}

// shapevec is a vector where each point has an angle greater than the last
std::complex<float> interceptAngle(float angle, std::vector<std::complex<float>> shapeVec) {
  for(int i = 0; i < shapeVec.size(); i++) {
    int lastIndex = i-1;
    if(i == 0) {
      lastIndex = shapeVec.size()-1;
    }
    std::complex<float> point2 = std::exp(complex_i * angle);
    std::complex<float> point3 = shapeVec[lastIndex];
    std::complex<float> point4 = shapeVec[i];
    std::pair<bool, float> intercept = interceptLines(std::complex<float>(0,0), point2, point3, point4);
    if(intercept.first && intercept.second <= 1 && intercept.second >= 0) {
      return intercept.second * point2;
    }
  }
  // this should not happen
  std::cout << "bad happened" << std::endl;
  return shapeVec[0];
  
}

std::vector<std::complex<float>> squareShapeVec = {std::complex<float>(-0.5, 0.5),
						   std::complex<float>(0.5, 0.5),
						   std::complex<float>(0.5, -0.5),
						   std::complex<float>(-0.5, -0.5)};
std::complex<float> squareFunction(float time) {
  
  time = std::fmod(time, 1);
  float angle = M_PI * 2 * time;
  return interceptAngle(angle, squareShapeVec);
}

struct Rotating {
  std::complex<float> coefficient;
  int n;
  std::complex<float> (*circleFunctionPointer)(float time);

  Rotating(int nIn, std::complex<float> coefficientIn) {
    n = nIn;
    coefficient = coefficientIn;
    circleFunctionPointer = circleFunction;
  }
};


struct State {
  bool mousePressed;
  std::vector<std::complex<float>> samples;
  std::vector<Rotating> rotatings;
  std::complex<float> (*circleFunctionPointer)(float time);
  

  State() {
    mousePressed = false;
    circleFunctionPointer = circleFunction;
  }
};

State state;


bool compareRotating(Rotating& r1, Rotating& r2) {
  return std::abs(r1.coefficient) > std::abs(r2.coefficient);
}



std::vector<Rotating> make_rotating(std::vector<std::complex<float>>& samples) {
  std::vector<Rotating> rotatings;
  for(int i = 0; i < granularity; i++) {
    int n = i + start;
    std::complex<float> c(0, 0);
    for(int samplei = 0; samplei < samples.size(); samplei++) {
      float t = ((float) samplei) / ((float) samples.size());
      c += samples[samplei] * state.circleFunctionPointer(-n * t);
    }
    c /= ((float) samples.size());
    Rotating r(n, c);
    r.circleFunctionPointer = state.circleFunctionPointer;
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


void renderLines() {
  

  glDisable(GL_LIGHTING);

  glColor3f(0, 0, 0);
  glLineWidth(3);
  glBegin(GL_LINES);
  
  float time = getTime() * 2 * M_PI * rotation_rate;
  std::complex<float> currentpos(0, 0);

  for(Rotating r : state.rotatings) {
    std::complex<float> oldpos(currentpos);
    std::complex<float> vec = r.circleFunctionPointer(r.n * time + std::arg(r.coefficient)/(2.0 * M_PI));
    vec *= std::abs(r.coefficient);
    currentpos += vec;

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

void printme(std::complex<float> a) {
  std::cout << a.real() << " and " << a.imag() << std::endl;
}

int main()
{
	initWindow();
	initGL();
	printHotKeys();
	renderLoop();
}
