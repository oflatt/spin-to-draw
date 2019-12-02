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
float fov = M_PI / 4.f;
float initialDistance = 4.5f;
float distance = 4.5f;

// per second
float rotation_rate = 0.05f;
float fovMoveSpeed = 0.1f;

struct Rotating {
  std::complex<float> coefficient;
  int n;
};

bool compareRotating(Rotating& r1, Rotating& r2) {
  return std::abs(r1.coefficient) > std::abs(r2.coefficient);
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

void computeNormals()
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
}

void loadObj(std::string p_path)
{
	std::ifstream nfile;
	nfile.open(p_path);
	std::string s;

	while (nfile >> s)
	{
		if (s.compare("v") == 0)
		{
			float x, y, z;
			nfile >> x >> y >> z;
			g_meshVertices.push_back(x);
			g_meshVertices.push_back(y);
			g_meshVertices.push_back(z);
		}
		else if (s.compare("f") == 0)
		{
			std::string sa, sb, sc;
			unsigned int a, b, c;
			nfile >> sa >> sb >> sc;

			a = std::stoi(sa);
			b = std::stoi(sb);
			c = std::stoi(sc);

			g_meshIndices.push_back(a - 1);
			g_meshIndices.push_back(b - 1);
			g_meshIndices.push_back(c - 1);
		}
		else
		{
			std::getline(nfile, s);
		}
	}

	computeNormals();

	std::cout << p_path << " loaded. Vertices: " << g_meshVertices.size() / 3 << " Triangles: " << g_meshIndices.size() / 3 << std::endl;
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

	float halfWidth = initialDistance * tan(fov / 2);

	// Perspective Projection
	if (enablePersp)
	{
		float currentFov = fov;

		distance = initialDistance;
		// Dolly zoom computation
		if (enableDolly) {
			currentFov = (fov - fmod((getTime() * M_PI * fovMoveSpeed), fov * 0.99));
			distance = halfWidth / tan(currentFov / 2);
		}


		float fovInDegree = radianToDegree(currentFov);
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

const std::complex<float> complex_i(0, 1);
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

void renderLines() {

  
  std::vector<std::complex<float>> samples = {std::complex<float>(-2.0, 2.0), std::complex<float>(0, 2.0),
						std::complex<float>(2.0, 2.0),
				std::complex<float>(2.0, 0.0), std::complex<float>(2.0, -2.0),
				std::complex<float>(0.0, -2.0), std::complex<float>(-2.0, -2.0)};

  std::vector<Rotating> rotatings = make_rotating(samples);
  
  

  glDisable(GL_LIGHTING);

  glColor3f(0, 0, 0);
  glLineWidth(3);
  glBegin(GL_LINES);
  
  float time = getTime() * 2 * M_PI * rotation_rate;
  std::complex<float> currentpos(0, 0);

  for(Rotating r : rotatings) {
    int n = r.n;
    std::complex<float> oldpos(currentpos);
    currentpos += r.coefficient * std::exp(n * 2 * M_PI * complex_i * time);
    glVertex3f(oldpos.real(), oldpos.imag(), distance);
    glVertex3f(currentpos.real(), currentpos.imag(), distance);
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
	renderLines();
	if (showCheckerboard)
		renderCheckerBoard();
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
	}
}

int main()
{
	initWindow();
	initGL();
	loadObj("data/teapot.obj");
	printHotKeys();
	renderLoop();
}
