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

int granularity = 50;
int start = -(granularity / 2);


unsigned int g_windowWidth = 1000;
unsigned int g_windowHeight = 800;
float modelWidth = 1.0;
float modelHeight = ((float) g_windowHeight) / ((float) g_windowWidth);
float lineZStep = 0.5 / granularity;
float number = 1;
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
bool zoomIn = false;
bool zoomOut = false;
bool goingUp = false;
bool goingDown = false;
bool goingLeft = false;
bool goingRight = false;

bool showCheckerboard = false;

// Dolly zoom options 
float fov = M_PI / 2.f;

// line up distance so that it matches the screen

float initialDistance = (1.0 / tan(fov/2)) * 0.25;
float distance = initialDistance;

// per second
float rotation_rate = 0.025f;

void renderDrawing();

const std::complex<float> complex_i(0, 1);

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
  float b = point2.real() - point1.real(); // 1.0
  float c = point4.real()- point3.real(); // 0.0
  float ap = point1.imag() - point3.imag();
  float bp = point2.imag()-point1.imag(); // 1.0
  float cp = point4.imag() - point3.imag(); // -1.0
  
  
  if (c == 0.0) {
	  if (b == 0.0) {
		  return std::make_pair(false, 0.0);
	  }
	  else {
		  return std::make_pair(true, -a / b);
	  }
  }
  
  if(c*bp - b*cp == 0.0) {
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
    std::complex<float> point1 = std::complex<float>(0, 0);
    std::complex<float> point2 = std::exp(complex_i * angle);
    std::complex<float> point3 = shapeVec[lastIndex];
    std::complex<float> point4 = shapeVec[i];
    //std::float mag = std::abs(point4 - point3);
    std::pair<bool, float> intercept = interceptLines(point3, point4, point1, point2);
    std::pair<bool, float> intercept2 = interceptLines(point1, point2, point3, point4);
    if(intercept.first && intercept.second <= 1.0  && intercept.second >= 0 && intercept2.second >= 0) {
      return intercept.second * (point4 - point3) + point3;
    }
  }
  // this should not happen
  std::cout << "bad happened" << std::endl;
  return shapeVec[0];
  
}

float squareSize = 1.0;
std::vector<std::complex<float>> squareShapeVec = {std::complex<float>(-squareSize, squareSize),
						   std::complex<float>(squareSize, squareSize),
						   std::complex<float>(squareSize, -squareSize),
						   std::complex<float>(-squareSize, -squareSize)};
std::complex<float> squareFunction(float time) {
  
  time = std::fmod(time, 1);
  float angle = -M_PI * 2 * time;
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
  double time = 0;
  float rotation;
  GLfloat g_modelViewMatrixState[16];
  std::complex<float> (*circleFunctionPointer)(float time);

  // state gets reset when click happens
  State() {
    rotation = 0;
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
      std::complex<float> circleInverse =  std::conj(state.circleFunctionPointer(-n * t));
      circleInverse *= 1.0/std::abs(circleInverse);
      c += samples[samplei] * circleInverse;
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

void zoom() {
	if (zoomOut) {
		float distance1 = -distance;
		distance = distance + 0.001;
		g_modelViewMatrix[14] = distance1;	
	}
	else if (zoomIn) {
		float distance1 = -distance;
		if (distance > 0.001) {
			distance = distance - 0.001;
		}
		g_modelViewMatrix[14] = distance1;
	}
	for (int i = 0; i < 16; i++) {
		state.g_modelViewMatrixState[i] = g_modelViewMatrix[i];
	}
	
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
  if(button == GLFW_MOUSE_BUTTON_LEFT) 
    {
      // reset state on button press
      if(action == GLFW_PRESS) {
	State newState;
	newState.circleFunctionPointer = state.circleFunctionPointer;
	state = newState;
	
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
void goingUpFunc() {
	if (goingUp) {
		number = number + 0.01;
		for (int i = 1; i < 16; i+=2) {
			if (i != 5) {
				g_modelViewMatrix[i] -= 0.001;
			}
		}
	}
	else if (goingDown) {
		number = number + 0.01;
		for (int i = 1; i < 16; i += 2) {
			if (i != 5) {
				g_modelViewMatrix[i] += 0.001;
			}
		}
	}
	else if (goingRight) {	
		number = number + 0.01;	
		for (int i = 0; i < 16; i += 2) {
			if (i != 5) {
				g_modelViewMatrix[i] -= 0.001;
			}
		}
	}
	else if (goingLeft) {
		number = number + 0.01;
		for (int i = 0; i < 16; i += 2) {
			if (i != 5) {
				g_modelViewMatrix[i] += 0.001;
			}
		}
	}
	for (int i = 0; i < 16; i++) {
		state.g_modelViewMatrixState[i] = g_modelViewMatrix[i];
	}
	
}
void glfwKeyCallback(GLFWwindow* p_window, int p_key, int p_scancode, int p_action, int p_mods)
{
  
  if (p_key == GLFW_KEY_ESCAPE && p_action == GLFW_PRESS)
    {
      glfwSetWindowShouldClose(g_window, GL_TRUE);
    }
  if (p_key == GLFW_KEY_RIGHT) {
    if(p_action == GLFW_PRESS) {
	  // Toggle Spinning
		if (teapotSpinLeft) {
			lastThingStopped = "left";
			teapotSpinLeft = false;
		}
		if (!teapotSpinRight) {
		  std::cout << "Spinning Right\n";
		  teapotSpinRight = true;
		}
	  
	  
    } else if(p_action == GLFW_RELEASE) {
      if(teapotSpinRight) {
		teapotSpinRight = false;
		lastThingStopped = "right";
      }
	  if (teapotSpinLeft) {
		  teapotSpinLeft = false;
		  lastThingStopped = "left";
	  }
    }
  }
  if (p_key == GLFW_KEY_UP) {
	  if (zoomIn) {
		  zoomIn = false;
	  }
	  else {
		  zoomIn = true;
	  }
	  if (p_action == GLFW_RELEASE) {
		  zoomIn = false;
	  }
	  zoomOut = false;
  }
  if (p_key == GLFW_KEY_DOWN) {
	  
	  if (zoomOut) {
		  zoomOut = false;
	  }
	  else {
		  zoomOut = true;
	  }
	  if (p_action == GLFW_RELEASE) {
		  zoomOut = false;
	  }
	  zoomIn = false;
  }
  if (p_key == GLFW_KEY_W) {
	  if (goingUp) {
		  goingUp = false;
	  }
	  else {
		  goingUp = true;
	  }
	  if (p_action == GLFW_RELEASE) {
		  goingUp = false;
	  }
	  goingDown = false;
  }
  if (p_key == GLFW_KEY_S) {
	  if (goingDown) {
		  goingDown = false;
	  }
	  else {
		  goingDown = true;
	  }
	  if (p_action == GLFW_RELEASE) {
		  goingDown = false;
	  }
	  goingUp = false;
  }
  if (p_key == GLFW_KEY_A) {
	  if (goingLeft) {
		  goingLeft = false;
	  }
	  else {
		  goingLeft = true;
	  }
	  if (p_action == GLFW_RELEASE) {
		  goingLeft = false;
	  }
	  goingRight = false;
  }
  if (p_key == GLFW_KEY_D) {
	  if (goingRight) {
		  goingRight = false;
	  }
	  else {
		  goingRight = true;
	  }
	  if (p_action == GLFW_RELEASE) {
		  goingRight = false;
	  }
	  goingLeft = false;
  }
  if (p_key == GLFW_KEY_LEFT) {
    if(p_action == GLFW_PRESS) {
	  // Toggle Spinning
		if (teapotSpinRight) {
			lastThingStopped = "right";
			teapotSpinRight = false;
		}
		if (!teapotSpinLeft) {
		  std::cout << "Spinning Left\n";
		  teapotSpinLeft = true;
		}
	  
    } else if(p_action == GLFW_RELEASE) {
      if(teapotSpinLeft) {
		teapotSpinLeft = false;
		lastThingStopped = "left";
		std::cout << "Teapot spinning off\n";
      }
	  if (teapotSpinRight) {
		  teapotSpinRight = false;
		  lastThingStopped = "right";
	  }
    }
  }

  if(p_key == GLFW_KEY_S && p_action == GLFW_PRESS) {
    if(state.circleFunctionPointer == circleFunction) {
      std::cout << "Switching to square shape";
      state.circleFunctionPointer = squareFunction;
    } else {
      std::cout << "Switching to circle shape";
      state.circleFunctionPointer = circleFunction;
    }
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

  //glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_DEPTH_TEST);
  glShadeModel(GL_SMOOTH);
}

void printHotKeys() {
	std::cout << "\nControls\n"
		  << "Click and drag to draw a shape\n"
		  << "Spin using left/right arrow keys\n"
		<< "Zoom using the up/down arrow keys\n"
		<< "Move the camera up, down, left or right with WASD Keys (it's reversed if you think about it as moving\n"
		<< "   the figure up, down, left or right)/n"
		<< "Change the shape of the vector rotation (circle or square fourier transform)\n"
		<< "Exit:                   Esc\n\n";
}

void clearModelViewMatrix()
{
	for (int i = 0; i < 4; ++i)
	{
		for (int j = 0; j < 4; ++j)
		{
			if ((4 * i + j) != 5.0) {
				g_modelViewMatrix[4 * i + j] = 0.0f;
			}
			
		}
	}
}

void updateModelViewMatrix()
{
	//clearModelViewMatrix();
	
	//state.time = abs(state.time - getTime());
	float rotation = state.rotation;
	if (lastThingStopped == "left") {

		g_modelViewMatrix[0] = cos(1 * rotation);
		g_modelViewMatrix[2] = sin(1 * rotation);
		//g_modelViewMatrix[5] = 1;
		g_modelViewMatrix[8] = sin(1 * rotation);
		g_modelViewMatrix[10] = -cos(1 * rotation);
	}
	else {
		g_modelViewMatrix[0] = cos(1 * rotation);
		g_modelViewMatrix[2] = -sin(1 * rotation);
		//g_modelViewMatrix[5] = 1;
		g_modelViewMatrix[8] = sin(1 * rotation);
		g_modelViewMatrix[10] = cos(1 * rotation);
	}

	if (teapotSpinLeft) {
		state.rotation = state.time * 2 * M_PI * rotation_rate;
		float rotation = state.rotation;
		g_modelViewMatrix[0] = cos(1 * rotation);
		g_modelViewMatrix[2] = sin(1 * rotation);
		//g_modelViewMatrix[5] = 1;
		g_modelViewMatrix[8] = sin(1 * rotation);
		g_modelViewMatrix[10] = -cos(1 * rotation);
		state.time = state.time + 0.25;
		//state.time = getTime()-state.time;
	}
	else if (teapotSpinRight) {
		state.rotation = (state.time * 2 * M_PI * rotation_rate);
		float rotation = state.rotation;
		g_modelViewMatrix[0] = cos(1 * rotation);
		g_modelViewMatrix[2] = -sin(1 * rotation);
		//g_modelViewMatrix[5] = 1;
		g_modelViewMatrix[8] = sin(1 * rotation);
		g_modelViewMatrix[10] = cos(1 * rotation);
		state.time = state.time  + 0.25;

	}

	g_modelViewMatrix[14] = -distance;
	g_modelViewMatrix[15] = 1.0f;
	//state.time = abs(state.time - getTime());

}
void setModelViewMatrix()
{
	glMatrixMode(GL_MODELVIEW);
	updateModelViewMatrix();
	zoom();
	goingUpFunc();
	glLoadMatrixf(g_modelViewMatrix);
}


void renderLines() {
  

  

  glColor3f(0, 0, 0);
  glLineWidth(3);
  glBegin(GL_LINES);
  
  float time = getTime() * 2 * M_PI * rotation_rate;
  std::complex<float> currentpos(0, 0);

  int i = 0;
  for(Rotating r : state.rotatings) {
    
    std::complex<float> oldpos(currentpos);
    std::complex<float> vec = r.coefficient * r.circleFunctionPointer(r.n * time);
    
    //vec *= std::abs(r.coefficient);
    currentpos += vec;

    std::complex<float> oldposScreen = model_to_screen(oldpos);
    std::complex<float> currentposScreen = model_to_screen(currentpos);
    glVertex3f(oldposScreen.real(), oldposScreen.imag(), -((int) state.rotatings.size()-i)*lineZStep);
    glVertex3f(currentposScreen.real(), currentposScreen.imag(), -((int) state.rotatings.size()-i)*lineZStep + lineZStep);
    i++;
  }

  glEnd();

}


void renderDrawing() {

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
	//g_modelViewMatrix[5] = 1;
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
  std::cout << interceptLines(std::complex<float>(0,0), std::complex<float>(-1.0,0.0),
			      std::complex<float>(0.5,1.0), std::complex<float>(0.5,-1.0)).second << std::endl;
	g_modelViewMatrix[5] = 1;
	initWindow();
	initGL();
	printHotKeys();
	renderLoop();
}
