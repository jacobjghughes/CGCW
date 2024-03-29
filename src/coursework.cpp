#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>
#include <glm/ext.hpp>
#include <CanvasPoint.h>
#include <CanvasTriangle.h>
#include <Colour.h> 
#include <TextureMap.h>
#include <TexturePoint.h>
#include <ModelTriangle.h>
#include <RayTriangleIntersection.h>
#include <fstream>
#include <string>
#include <iostream>
#include <unordered_map>
#include <cmath>
#include <algorithm>

#define WIDTH 600 //300
#define HEIGHT 600 //300
// #define BRICK "../textures/texture.ppm"
#define MODELS "./models"
#define BOX "./models/textured-cornell-box.obj"
#define BALL "./models/sphere.obj"
#define LOGO "./models/logo1.obj"
// #define BOX "/home/jakehughes/Desktop/CG2021/Weekly Workbooks/04 Wireframes and Rasterising/Wireframes/models/boxes.obj"
#define MATS "./models/textured-cornell-box.mtl"
#define FOCLEN 6 // default focal length
#define SCALE 50
#define DEFPOS glm::vec3(0.0, 0.0, 1.0)
// #define DEFPOS glm::vec3(0.0, 0.0, 10.0)
#define DEFROT glm::mat3(1,0,0,0,1,0,0,0,1)
#define ORIGIN CanvasPoint(0,0,0)
// #define DEFLIGHT {glm::vec3(0.4, 0.4, 0)}
#define DEFINT 0.7
#define DEFCOL Colour(255,255,255)
#define DEFCOLINT 0.8
// #define DEFLIGHT {Light(glm::vec3(0, 0.4, 0), 0.5, Colour(255,0,255), 0.2), Light(glm::vec3(-0.2, -0.3, 0.5), 0.7, Colour(255, 200, 200), 0.2)}
#define DEFLIGHT {Light(glm::vec3(0, 0.4, 0), 0.7, Colour(255, 200, 200), 0.1), Light(glm::vec3(-0.2, -0.3, 0.5), 0.5, Colour(128,30,180), 0.3)}
#define SPECINT 200
#define SPECPOW 1024
#define PROXINT 6

// camerawork stuff ------------------------------------------------------------------------------
	std::vector<glm::vec3> positions({
		(glm::vec3(0, 0, 1)),
		(glm::vec3(-0.282683, -0.2, 0.623879)),
		(glm::vec3(-0.0224172, 0.408267, 0.621784)),
		(glm::vec3(0.4, 0.294236, -0.0585273)),
		(glm::vec3(0.4, 0.294236, -0.0585273))
	});

	std::vector<glm::mat3> rotations({
		(glm::mat3(1, 0, 0, 0, 1, 0, 0, 0, 1)),
		(glm::mat3(0.852415, -0, 0.522866, 0, 1, 0, -0.522866, 0, 0.852415)),
		(glm::mat3(0.999351, -0, 0.0360296, 0.0227939, 0.774443, -0.632234, -0.0279028, 0.632644, 0.77394)),
		(glm::mat3(-0.144777, 0, -0.989464, -0.582272, 0.808518, 0.085197, 0.8, 0.588472, -0.117055)),
		(glm::mat3(-0.144777, 0, -0.989464, -0.582272, 0.808518, 0.085197, 0.8, 0.588472, -0.117055))
	});

	std::vector<Colour> L1Col({Colour(109, 18, 170), Colour(146, 237, 85)});
	std::vector<Colour> L2Col({Colour(115, 13, 14), Colour(140, 242, 241)});
	// std::vector<float> L1Int({0.1, 0.5, 0.8, 0.9, 0.9});
	// std::vector<float> L2Int({0.8, 0.8, 0.2, 0.0, 0.0});
	std::vector<float> L1Int({0.0, 0.0, 1.0, 1.0, 0.0, 0.8});
	std::vector<float> L2Int({0.0, 1.0, 1.0, 0.0, 0.0, 0.8});
	// std::vector<int> keyframeTime({3, 3, 3, 1});
	// std::vector<int> keyframeTime({30, 30, 30, 15});
	std::vector<int> keyframeTime({20, 20, 20, 20, 20, 60});
	int keyFrame = 0;
	int frame = 0;
	glm::vec3 posstep;
	glm::mat3 rotstep;
	float l1step;
	float l2step;
	float l1colstepr, l1colstepg, l1colstepb;
	float l2colstepr, l2colstepg, l2colstepb;

// -------------------------------------------------------------------------

struct Light {
	glm::vec3 light;
	float intensity{DEFINT};
	Colour colour{DEFCOL};
	float colourIntensity{DEFCOLINT};


	Light() {
		light = glm::vec3(0,0,0);
	}

	Light(glm::vec3 vec, float i = DEFINT, Colour c = DEFCOL, float ci = DEFCOLINT) {
		light = vec;
		intensity = i;
		colour = c;
		colourIntensity = ci;
	}
};

bool OrbiterToggle;
glm::vec3 campos(DEFPOS); 
float depthBuffer[WIDTH * HEIGHT];
std::vector<Light> lightList = DEFLIGHT;
std::vector<ModelTriangle> triangleList;
int modelCount = 0;
std::vector<TextureMap> textureList;
glm::mat3 camrot(DEFROT);
int renderStyle = 4;

CanvasPoint vecToCP(glm::vec3 vec) {
	CanvasPoint p(vec.x, vec.y, vec.z);
	return p;
}

glm::vec3 cpToVec(CanvasPoint p) {
	glm::vec3 vec(p.x, p.y, p.depth);
	return vec;
}

CanvasTriangle modelToCanvas(ModelTriangle mt) {
	CanvasTriangle ct;
	for (int i = 0; i < 3; i++) {
		ct.vertices[i] = vecToCP(mt.vertices[i]);
		ct.vertices[i].texturePoint = mt.texturePoints[i];
	}
	return ct;
}

Colour scaleColour(float scale, Colour c) {
	Colour nc;
	nc.red = c.red * scale;
	nc.blue = c.blue * scale;
	nc.green = c.green * scale;
	return nc;
}

glm::vec3 cpToLowPVec(CanvasPoint p) {
	glm::lowp_vec3 vec(p.x, p.y, p.depth);
	return vec;
}

uint32_t ColourToInt(Colour c) {
	uint32_t colour = (255 << 24) + (int(c.red) << 16) + (int(c.green) << 8) + int(c.blue);
	return colour;
}

bool inbounds(CanvasPoint p) {
	if (p.x <= 0 || p.x >= WIDTH) return false;
	else if (p.y <= 0 || p.y >= HEIGHT) return false;
	else if (p.depth < 0) return false; 
	else return true;
}

bool inbounds(int x, int y) {
	if (x <= 0 || x >= WIDTH) return false;
	else if (y <= 0 || y >= HEIGHT) return false;
	else return true;
}

bool inrange(int x, int y, float min, float max) {
	if (x < min || x > max) return false;
	else if (y < min || y > max) return false;
	else return true;
}

bool verifyTexturePoints(std::array<TexturePoint, 3> points) {
	for (int i = 0; i < 3; i++) {
		if ((points[i].x > 0 || points[i].y > 0)) {
			std::cout << "index: " << i << " tpx: " << points[i].x << " tpy: " << points[i].y << std::endl;
			return true;
		}
	}
	return false;
}

float clamp(float x, float min, float max) {
	if (x <= min) return min;
	else if (x >= max) return max;
	else return x;
}

std::vector<glm::lowp_vec3> interpolateLowPVec(glm::lowp_vec3 from, glm::lowp_vec3 to, int numberOfValues) {
	std::vector<glm::lowp_vec3> interpolatedList;
	glm::lowp_vec3 step = (to - from);

	if (numberOfValues < 2) numberOfValues = 2;
	step /= (numberOfValues - 1);

	interpolatedList.push_back(from);
	for (int i = 1; i < numberOfValues; i++) {
		if (interpolatedList.back().x >= WIDTH || interpolatedList.back().x <= 0 || interpolatedList.back().y <= 0 || interpolatedList.back().y >= HEIGHT) {
			std::cout << "x: " << interpolatedList.back().x << " y: " << interpolatedList.back().y << " z: " << interpolatedList.back().z << std::endl;
			std::cout << "from: " << vecToCP(from) << " to: " << vecToCP(to) << std::endl;
			i = numberOfValues;
		}
		interpolatedList.push_back(interpolatedList.back() + step);
	}
	return interpolatedList;
}

std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, int numberOfValues) {
	std::vector<glm::vec3> interpolatedList;
	glm::vec3 step = (to - from);

	if (numberOfValues < 2) numberOfValues = 2;
	step /= (numberOfValues - 1);

	interpolatedList.push_back(from);
	for (int i = 1; i < numberOfValues; i++) {
		interpolatedList.push_back(interpolatedList.back() + step);
	}
	return interpolatedList;
}

std::vector<CanvasPoint> interpolateTexturePoints(CanvasPoint from, CanvasPoint to, int numberOfValues) {
	std::vector<CanvasPoint> interpolatedList;
	float xstep, ystep, zstep, txstep, tystep;

	if (numberOfValues < 2) numberOfValues = 2;
	xstep = (to.x - from.x) / (numberOfValues - 1);
	ystep = (to.y - from.y) / (numberOfValues - 1);
	zstep = (to.depth - from.depth) / (numberOfValues - 1);
	txstep = (to.texturePoint.x - from.texturePoint.x) / (numberOfValues - 1);
	tystep = (to.texturePoint.y - from.texturePoint.y) / (numberOfValues - 1);

	interpolatedList.push_back(from);
	for (int i = 1; i < numberOfValues; i++) {
		CanvasPoint p(interpolatedList.back());
		p.x += xstep;
		p.y += ystep;
		p.depth += zstep;
		p.texturePoint.x += txstep;
		p.texturePoint.y += tystep;
		interpolatedList.push_back(p);
	}
	return interpolatedList;
}

void drawLine(CanvasPoint from, CanvasPoint to, DrawingWindow &window, Colour c = Colour(230,10,230)) {
	int xDiff = round(abs(to.x - from.x)); 
	int yDiff = round(abs(to.y - from.y));
	int numberOfValues = std::max(xDiff, yDiff);
	numberOfValues *= 2;
	uint32_t colour = ColourToInt(c);

	glm::vec3 vecfrom = cpToVec(from);
	glm::vec3 vecto = cpToVec(to);

	std::vector<glm::vec3> Points = interpolateThreeElementValues(vecfrom, vecto, numberOfValues);
	for (int i = 0; i < Points.size(); i++) {
		int x = round(Points[i].x);
		int y = round(Points[i].y);
		float z = Points[i].z;
		int index = y * WIDTH + x;

		// only draw pixel if its on canvas, and checks index is within range
		if (inbounds(CanvasPoint(x,y)) && index < WIDTH*HEIGHT && index >= 0) {
			if (z > 0 && (depthBuffer[index] < 100/z || depthBuffer[index] == 0)) {
				depthBuffer[index] = 100/z;
				window.setPixelColour(x, y, colour);
			}	
		}
	}
}

void drawStrokedTriangle(ModelTriangle t, DrawingWindow &window, Colour c = Colour(230,10,230)) {
	drawLine(vecToCP(t.vertices[0]), vecToCP(t.vertices[1]), window, c);
	drawLine(vecToCP(t.vertices[0]), vecToCP(t.vertices[2]), window, c);
	drawLine(vecToCP(t.vertices[1]), vecToCP(t.vertices[2]), window, c);
	
	// drawLine(t.v0(), t.v1(), window, c);
	// drawLine(t.v0(), t.v2(), window, c);
	// drawLine(t.v1(), t.v2(), window, c);
}

void drawStrokedTriangle(CanvasTriangle t, DrawingWindow &window, Colour c = Colour(230,10,230)) {
	drawLine(t.v0(), t.v1(), window, c);
	drawLine(t.v0(), t.v2(), window, c);
	drawLine(t.v1(), t.v2(), window, c);
}

CanvasPoint getCanvasIntersectionPoint(glm::vec3 vertexPosition, float focalLength = FOCLEN, float scaleFactor = SCALE) {
	// transpose vertex position so model vertex is in camera coordinate system
	vertexPosition -= campos;
	vertexPosition = vertexPosition * camrot;

	CanvasPoint projectedPoint;
	projectedPoint.x = -1 * SCALE * (focalLength * vertexPosition.x / vertexPosition.z);
	projectedPoint.y = SCALE * (focalLength * vertexPosition.y / vertexPosition.z);
	// projectedPoint = scalePoint(projectedPoint, scaleFactor);
	// projectedPoint.x = WIDTH - projectedPoint.x;
	projectedPoint.x += (WIDTH / 2);
	projectedPoint.y += (HEIGHT / 2);
	projectedPoint.depth = -1 * focalLength * vertexPosition.z;
	// std::cout << "depth: " << projectedPoint.depth << std::endl;
	return projectedPoint;
}

ModelTriangle sortVerticesVertically(ModelTriangle t) {
	// add triangle points to list of vertices
	std::vector<glm::vec3> vertices;
	vertices.push_back(t.vertices[0]);
	vertices.push_back(t.vertices[1]);
	vertices.push_back(t.vertices[2]);

	// sort vertices vertically 
	if (t.vertices[0].y > t.vertices[1].y) {
		std::swap(t.vertices[0],t.vertices[1]);
		std::swap(t.texturePoints[0], t.texturePoints[1]);
		
		if (t.vertices[0].y > t.vertices[2].y) {
			std::swap(t.vertices[0],t.vertices[2]);
			std::swap(t.texturePoints[0], t.texturePoints[2]);
		}
		if (t.vertices[1].y > t.vertices[2].y) {
			std::swap(t.vertices[1],t.vertices[2]);
			std::swap(t.texturePoints[1], t.texturePoints[2]);
		}
	}

	else if (t.vertices[0].y > t.vertices[2].y) {
		std::swap(t.vertices[0],t.vertices[2]);
		std::swap(t.texturePoints[0], t.texturePoints[2]);

		if (t.vertices[1].y > t.vertices[2].y) {
			std::swap(t.vertices[1],t.vertices[2]);
			std::swap(t.texturePoints[1], t.texturePoints[2]);
		}
	}

	else if (t.vertices[1].y > t.vertices[2].y) {
		std::swap(t.vertices[1],t.vertices[2]);
		std::swap(t.texturePoints[1], t.texturePoints[2]);
	}
	return t;
}

CanvasTriangle sortVerticesVertically(CanvasTriangle t) {
	// add triangle points to list of vertices
	std::vector<CanvasPoint> vertices;
	vertices.push_back(t.v0());
	vertices.push_back(t.v1());
	vertices.push_back(t.v2());

	// sort vertices vertically 
	if (vertices[0].y > vertices[1].y) {
		std::swap(vertices[0],vertices[1]);
		
		if (vertices[0].y > vertices[2].y) {
			std::swap(vertices[0],vertices[2]);
		}
		if (vertices[1].y > vertices[2].y) {
			std::swap(vertices[1],vertices[2]);
		}
	}

	else if (vertices[0].y > vertices[2].y) {
		std::swap(vertices[0],vertices[2]);

		if (vertices[1].y > vertices[2].y) {
			std::swap(vertices[1],vertices[2]);
		}
	}

	else if (vertices[1].y > vertices[2].y) {
		std::swap(vertices[1],vertices[2]);
	}

	t = CanvasTriangle(vertices[0],vertices[1],vertices[2]);
	return t;
}

CanvasTriangle sortTopTriangle(CanvasTriangle t) {
	t = sortVerticesVertically(t);
	if(t.v1().x > t.v2().x) std::swap(t.v1(),t.v2());
	return t;
}

ModelTriangle sortTopTriangle(ModelTriangle t) {
	t = sortVerticesVertically(t);
	if(t.vertices[1].x > t.vertices[2].x) {
		std::swap(t.vertices[1],t.vertices[2]);
		std::swap(t.texturePoints[1], t.texturePoints[2]);
	}
	return t;
}

CanvasTriangle sortBottomTriangle(CanvasTriangle t) {
	t = sortVerticesVertically(t);
	std::swap(t.v2(), t.v0());
	if(t.v1().x > t.v2().x) std::swap(t.v1(),t.v2());
	return t;
}

ModelTriangle sortBottomTriangle(ModelTriangle t) {
	t = sortVerticesVertically(t);
	std::swap(t.vertices[2], t.vertices[0]);
	std::swap(t.texturePoints[2], t.texturePoints[0]);

	if(t.vertices[1].x > t.vertices[2].x) {
		std::swap(t.vertices[1],t.vertices[2]);
		std::swap(t.texturePoints[1], t.texturePoints[2]);
	}
	return t;
}

CanvasPoint lineRatio(CanvasPoint from, CanvasPoint to, float ratio) {
	float xlen = (to.x - from.x) * ratio;
	float ylen = (to.y - from.y) * ratio;
	float zlen = (to.depth - from.depth) * ratio;

	CanvasPoint newp(from);
	newp.x += xlen;
	newp.y += ylen;
	newp.depth += zlen;
	return newp;
}

TexturePoint lineRatio(TexturePoint from, TexturePoint to, float ratio) {
	float xlen = (to.x - from.x) * ratio;
	float ylen = (to.y - from.y) * ratio;

	TexturePoint newp(from);
	newp.x += xlen;
	newp.y += ylen;
	return newp;
}

CanvasPoint findIntercept(CanvasTriangle t) {
	float verticalHeight = (abs(t.v2().y - t.v0().y));
	float interceptHeight = (abs(t.v1().y - t.v0().y));
	float ratio = interceptHeight / verticalHeight;

	CanvasPoint intercept(lineRatio(t.v0(), t.v2(), ratio));
	return intercept;
}

CanvasPoint findIntercept(ModelTriangle t) {
	float verticalHeight = (abs(t.vertices[2].y - t.vertices[0].y));
	float interceptHeight = (abs(t.vertices[1].y - t.vertices[0].y));
	float ratio = interceptHeight / verticalHeight;

	CanvasPoint intercept(lineRatio(vecToCP(t.vertices[0]), vecToCP(t.vertices[2]), ratio));
	intercept.texturePoint = TexturePoint(lineRatio(t.texturePoints[0], t.texturePoints[2], ratio));
	return intercept;
}

TexturePoint findTexturePointIntercept(CanvasTriangle t) {
	float verticalHeight = (abs(t.v2().y - t.v0().y));
	float interceptHeight = (abs(t.v1().y - t.v0().y));
	float ratio = interceptHeight / verticalHeight;

	TexturePoint intercept(lineRatio(t.v0().texturePoint, t.v2().texturePoint, ratio));
	return intercept;
}

TexturePoint findTexturePointIntercept(ModelTriangle t) {
	float verticalHeight = (abs(t.vertices[2].y - t.vertices[0].y));
	float interceptHeight = (abs(t.vertices[1].y - t.vertices[0].y));
	float ratio = interceptHeight / verticalHeight;

	TexturePoint intercept(lineRatio(t.texturePoints[0], t.texturePoints[2], ratio));
	return intercept;
}

void rasterizeTriangle(CanvasTriangle t, DrawingWindow &window, Colour c) {

	t = sortVerticesVertically(t); 
	// drawStrokedTriangle(t, window, c);

	std::vector<CanvasPoint> vertices;
	vertices.push_back(t.v0());
	vertices.push_back(t.v1());
	vertices.push_back(t.v2());


	CanvasTriangle topT; 
	CanvasTriangle botT;
	
	bool flatTop = false;
	bool flatBot = false;
	// checks the triangle to see if it already has a flat top or bottom
	// bottom triangle (flat top)
	if (floor(vertices[0].y) == floor(vertices[1].y)) {
		flatTop = true;
		botT = sortBottomTriangle(t);
	}

	// top triangle (flat bottom)
	if (floor(vertices[1].y) == floor(vertices[2].y)) {
		flatBot = true;
		topT = sortTopTriangle(t);
	}

	if (!flatTop && !flatBot) {
		// creates new point intercept which is the horizontal intercept at the y value of the middle vertex
		CanvasPoint intercept(findIntercept(t));
		botT = sortBottomTriangle(CanvasTriangle(intercept, vertices[1], vertices[2]));
		topT = sortTopTriangle(CanvasTriangle(intercept, vertices[0], vertices[1]));
	}

	//rasterise top triangle, if there is one
	if (!flatTop) {

		int height = abs(round(topT.v0().y) - round(topT.v1().y)) * 2;
		int base = abs(round(topT.v1().x) - round(topT.v2().x));
		int left = abs(round(topT.v1().x) - round(topT.v0().x));
		int right = abs(round(topT.v0().x - round(topT.v2().x)));

		int width = std::max(left, right);
		width = std::max(width, base);
	
		std::vector<glm::vec3>  ApexLeftY = interpolateThreeElementValues(cpToVec(topT.v0()), cpToVec(topT.v1()), height);
		std::vector<glm::vec3> ApexRightY = interpolateThreeElementValues(cpToVec(topT.v0()), cpToVec(topT.v2()), height);

		std::vector<glm::vec3>  LeftApexX = interpolateThreeElementValues(cpToVec(topT.v1()), cpToVec(topT.v0()), left);
		std::vector<glm::vec3> ApexRightX = interpolateThreeElementValues(cpToVec(topT.v0()), cpToVec(topT.v2()), right);
		std::vector<glm::vec3> LeftRight = interpolateThreeElementValues(cpToVec(topT.v1()), cpToVec(topT.v2()), base);
		
		// horizontal rake row
		for (int ys = 0; ys < height; ys++) {
			CanvasPoint from, to;
			from = vecToCP(ApexLeftY[ys]);
			to  = vecToCP(ApexRightY[ys]);
			drawLine(from, to, window, c);
		}
		
		// vertical rake row
		for (int xs = 0; xs < width; xs++) {
			CanvasPoint from, to;

			// if top triangle is obtuse draw vertical line between [v0,v1]@xs and [v0,v2]@xs
			if (width > base) {
				// if top triangle is left leaning
				if (round(topT.v0().x) < round(topT.v1().x)) {

					// if xs is to the left of base
					if (xs < left) {
						from = vecToCP(LeftApexX[left - xs - 1]);
						to = vecToCP(ApexRightX[xs]);
					}

					// else xs is above the base
					else {
						from = vecToCP(LeftRight[xs - left]);
						to = vecToCP(ApexRightX[xs]);
					}
				}
				// else if top triangle is right leaning
				else if (round(topT.v0().x) > round(topT.v2().x)) {
					
					// if xs is above the base
					if (xs < base) {
						from = vecToCP(LeftApexX[xs]);
						to = vecToCP(LeftRight[xs]);
					}
					// else xs is to the right of the base
					else {
						from = vecToCP(LeftApexX[xs]);
						to = vecToCP(ApexRightX[right - (xs - base) - 1]);
					}

				}
			}
			
			// else top triangle is acute
			else {
				if (xs < left) {
					from = vecToCP(LeftApexX[xs]);
					to = vecToCP(LeftRight[xs]);
				}
				else {
					from = vecToCP(ApexRightX[xs - left]);
					to = vecToCP(LeftRight[xs]);
				}
			}

			drawLine(from, to, window, c);
		}
	}
	// rasterise bottom triangle, if there is one
	if (!flatBot) {
		int height = abs(round(botT.v0().y) - round(botT.v1().y)) * 2;
		int base = abs(round(botT.v1().x) - round(botT.v2().x));
		int left = abs(round(botT.v1().x) - round(botT.v0().x));
		int right = abs(round(botT.v0().x - round(botT.v2().x)));

		int width = std::max(left, right);
		width = std::max(width, base);

		std::vector<glm::vec3>  ApexLeftY = interpolateThreeElementValues(cpToVec(botT.v0()), cpToVec(botT.v1()), height);
		std::vector<glm::vec3> ApexRightY = interpolateThreeElementValues(cpToVec(botT.v0()), cpToVec(botT.v2()), height);

		std::vector<glm::vec3>  LeftApexX = interpolateThreeElementValues(cpToVec(botT.v1()), cpToVec(botT.v0()), left);
		std::vector<glm::vec3> ApexRightX = interpolateThreeElementValues(cpToVec(botT.v0()), cpToVec(botT.v2()), right);
		std::vector<glm::vec3> LeftRight = interpolateThreeElementValues(cpToVec(botT.v1()), cpToVec(botT.v2()), base);
	
		// horizontal rake row
		for (int ys = 0; ys < height; ys++) {
			CanvasPoint from, to;
			from = vecToCP(ApexLeftY[ys]);
			to  = vecToCP(ApexRightY[ys]);
			drawLine(from, to, window, c);
		}
		
		// vertical rake row
		for (int xs = 0; xs < width; xs++) {
			CanvasPoint from, to;

			// if bot triangle is obtuse draw vertical line between [v0,v1]@xs and [v0,v2]@xs
			if (width > base) {
				// if bot triangle is left leaning
				if (round(botT.v0().x) < round(botT.v1().x)) {

					// if xs is to the left of base
					if (xs < left) {
						from = vecToCP(LeftApexX[left - xs - 1]); // mod left?
						to = vecToCP(ApexRightX[xs]);
					}

					// else xs is above the base
					else {
						from = vecToCP(LeftRight[xs - left]);
						to = vecToCP(ApexRightX[xs]);
					}
				}
				// else if top triangle is right leaning
				else if (round(botT.v0().x) > round(botT.v2().x)) {
					
					// if xs is above the base
					if (xs < base) {
						from = vecToCP(LeftApexX[xs]);
						to = vecToCP(LeftRight[xs]);
					}
					// else xs is to the right of the base
					else {
						from = vecToCP(LeftApexX[xs]);
						to = vecToCP(ApexRightX[right - (xs - base) - 1]); // mod right??
					}

				}
			}
			
			// else top triangle is acute
			else {
				if (xs < left) {
					from = vecToCP(LeftApexX[xs]);
					to = vecToCP(LeftRight[xs]);
				}
				else {
					from = vecToCP(ApexRightX[xs - left]);
					to = vecToCP(LeftRight[xs]);
				}
			}

			drawLine(from, to, window, c);
		}
	}
}

uint32_t getTexturePixel(float x, float y, int textureIndex) {
	TextureMap map = textureList.at(textureIndex);
	return map.pixels[(round(y) * map.width) + round(x)];
}

uint32_t getTexturePixel(TexturePoint p, TextureMap &map) {
	return map.pixels[(round(p.y) * map.width) + round(p.x)];
}

void textureTriangle(ModelTriangle t, DrawingWindow &window, int textureIndex = 0) {
	TextureMap texture = textureList.at(textureIndex);
	// convert texture percentages to texture coordinates from loaded texture 
	for (int i = 0; i < 3; i++) {
		t.texturePoints[i].x *= texture.width;
		t.texturePoints[i].y *= texture.height;
	}
	
	// std::cout << texture << std::endl;
	// MAKE A TEXTUREMAP DICTIONARY

	t = sortVerticesVertically(t); 

	// std::vector<CanvasPoint> vertices;
	// vertices.push_back(t.v0());
	// vertices.push_back(t.v1());
	// vertices.push_back(t.v2());

	ModelTriangle topT; 
	ModelTriangle botT;

	bool flatTop = false;
	bool flatBot = false;

	// checks the triangle to see if it already has a flat top or bottom
	// bottom triangle (flat top)
	if (floor(t.vertices[0].y) == floor(t.vertices[1].y)) {
		flatTop = true;
		botT = sortBottomTriangle(t);
	}

	// top triangle (flat bottom)
	if (floor(t.vertices[1].y) == floor(t.vertices[2].y)) {
		flatBot = true;
		topT = sortTopTriangle(t);
	}

	if (!flatTop && !flatBot) {
		// creates new point intercept which is the horizontal intercept at the y value of the middle vertex
		CanvasPoint intercept(findIntercept(t));
		// intercept.texturePoint = findTexturePointIntercept(t);

		topT = ModelTriangle(cpToVec(intercept), t.vertices[0], t.vertices[1], t.colour);
		topT.texturePoints = {{intercept.texturePoint, t.texturePoints[0], t.texturePoints[1]}};
		topT = sortTopTriangle(topT);

		botT = ModelTriangle(cpToVec(intercept), t.vertices[1], t.vertices[2], t.colour);
		botT.texturePoints = {{intercept.texturePoint, t.texturePoints[1], t.texturePoints[2]}};
		botT = sortBottomTriangle(botT);
	}

	if (!flatTop) {	
		// Interpolate left and right points on top triangle and corresponding texturepoints
		int height = abs(round(topT.vertices[1].y - topT.vertices[0].y));// * 2;
		CanvasPoint Apex, Left, Right;
		Apex = vecToCP(topT.vertices[0]);
		Apex.texturePoint = topT.texturePoints[0];

		Left = vecToCP(topT.vertices[1]);
		Left.texturePoint = topT.texturePoints[1];

		Right = vecToCP(topT.vertices[2]);
		Right.texturePoint = topT.texturePoints[2];

		std::vector<CanvasPoint> ApexLeftY  = interpolateTexturePoints(Apex, Left, height);
		std::vector<CanvasPoint> ApexRightY = interpolateTexturePoints(Apex, Right, height);

		// draw horizontal lines between [apex,left]@ys and [apex,right]@ys
		for (int ys = 0; ys < height; ys++) {
			CanvasPoint from, to; 
			int width = abs(round(ApexLeftY[ys].x - ApexRightY[ys].x));
			// interpolate rake row on texture points
			std::vector<CanvasPoint> rakeRow = interpolateTexturePoints(ApexLeftY[ys],ApexRightY[ys], width);
			for (int xs = 0; xs < width; xs++) {
				CanvasPoint p;
				p.x = round(rakeRow[xs].x);
				p.y = round(rakeRow[xs].y);
				p.depth = round(rakeRow[xs].depth);
				p.texturePoint.x = round(rakeRow[xs].texturePoint.x);
				p.texturePoint.y = round(rakeRow[xs].texturePoint.y);
				int index = p.y * WIDTH + p.x;

				if (inbounds(p) && index < WIDTH*HEIGHT && index >= 0) {
					if (p.depth > 0 && (depthBuffer[index] < 100 / p.depth || depthBuffer[index] == 0)) {
						depthBuffer[index] = 100/p.depth;
						window.setPixelColour(round(p.x), round(p.y), getTexturePixel(p.texturePoint, texture));
					}	
				}
			}
		}
	}

	if (!flatBot) {
		// Interpolate left and right points on bot triangle and corresponding texturepoints
		int height = abs(round(botT.vertices[1].y - botT.vertices[0].y));// * 2;
		CanvasPoint Apex, Left, Right;
		Apex = vecToCP(botT.vertices[0]);
		Apex.texturePoint = botT.texturePoints[0];

		Left = vecToCP(botT.vertices[1]);
		Left.texturePoint = botT.texturePoints[1];

		Right = vecToCP(botT.vertices[2]);
		Right.texturePoint = botT.texturePoints[2];

		std::vector<CanvasPoint> ApexLeftY  = interpolateTexturePoints(Apex, Left, height);
		std::vector<CanvasPoint> ApexRightY = interpolateTexturePoints(Apex, Right, height);


		// draw horizontal lines between [apex,left]@ys and [apex,right]@ys
		for (int ys = 0; ys < height; ys++) {
			CanvasPoint from, to; 
			int width = abs(round(ApexLeftY[ys].x - ApexRightY[ys].x));
			// interpolate rake row on texture points
			std::vector<CanvasPoint> rakeRow = interpolateTexturePoints(ApexLeftY[ys],ApexRightY[ys], width);
			for (int xs = 0; xs < width; xs ++) {
				CanvasPoint p;
				p.x = round(rakeRow[xs].x);
				p.y = round(rakeRow[xs].y);
				p.depth = round(rakeRow[xs].depth);
				p.texturePoint.x = round(rakeRow[xs].texturePoint.x);
				p.texturePoint.y = round(rakeRow[xs].texturePoint.y);
				int index = p.y * WIDTH + p.x;

				if (inbounds(p) && index < WIDTH*HEIGHT && index >= 0) {
					if (p.depth > 0 && (depthBuffer[index] < 100 / p.depth || depthBuffer[index] == 0)) {
						depthBuffer[index] = 100/p.depth;
						window.setPixelColour(round(p.x), round(p.y), getTexturePixel(p.texturePoint, texture));
					}	
				}
			}
		}
	}
	drawStrokedTriangle(t,window,Colour(255,255,255));
}


RayTriangleIntersection getClosestIntersection(glm::vec3 directionVector, glm::vec3 startpos = campos, int ignoreIndex = -1) {
	RayTriangleIntersection intersection;
	intersection.distanceFromCamera = INFINITY;

	for (int i = 0; i < triangleList.size(); i++) {
		
		if (i == ignoreIndex) continue;

		glm::vec3 e0 = triangleList[i].vertices[1] - triangleList[i].vertices[0];
		glm::vec3 e1 = triangleList[i].vertices[2] - triangleList[i].vertices[0];
		glm::vec3 SPVector = startpos - triangleList[i].vertices[0];
		glm::mat3 DEMatrix(-directionVector, e0, e1);
		glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;

		if (possibleSolution.x < intersection.distanceFromCamera && possibleSolution.x > 0 
			&& possibleSolution.y <= 1 && possibleSolution.y >= 0 
			&& possibleSolution.z <= 1 && possibleSolution.z >= 0 
			&& (possibleSolution.y + possibleSolution.z <= 1)) {

				intersection.distanceFromCamera = possibleSolution.x;
				intersection.intersectedTriangle = triangleList[i];
				intersection.intersectionPoint = glm::vec3(triangleList[i].vertices[0] + possibleSolution.y * e0 + possibleSolution.z * e1);
				intersection.triangleIndex = i;
				intersection.solution = possibleSolution;
		}
	}
	// if mirrored, then reflect again
	if (intersection.intersectedTriangle.reflectiveness == 1) {
		glm::vec3 reflectedRay;
		if (intersection.intersectedTriangle.shaded) {
			float u, v, w;
			u = intersection.solution.y;
			v = intersection.solution.z;
			w = 1 - u - v;
			glm::vec3 pointNormal;
			pointNormal += intersection.intersectedTriangle.vertexNormals[0] * w;
			pointNormal += intersection.intersectedTriangle.vertexNormals[1] * u;
			pointNormal += intersection.intersectedTriangle.vertexNormals[2] * v;
			pointNormal = glm::normalize(pointNormal);

			float angleDiff = glm::dot(directionVector, pointNormal);
			if (angleDiff < 0) {
				pointNormal = -pointNormal;
				angleDiff = glm::dot(directionVector, pointNormal);
			}

			reflectedRay = directionVector - 2 * pointNormal * (glm::dot(directionVector, pointNormal));
		}
		else {
			reflectedRay = directionVector - 2 * intersection.intersectedTriangle.normal * (glm::dot(directionVector, intersection.intersectedTriangle.normal));
		}
		intersection = getClosestIntersection(glm::normalize(reflectedRay), intersection.intersectionPoint, intersection.triangleIndex);
	}
	return intersection;
}

float specularLighting(RayTriangleIntersection &projectedPoint, glm::vec3 view, Light &thisLight, int reflectiveness = SPECPOW) {
	// cast ray from selected light to point
	glm::vec3 incidentRay = glm::normalize(projectedPoint.intersectionPoint - thisLight.light);
	// reflected ray is derived from following formula (Ri - 2 N (Ri . N))
	glm::vec3 reflectedRay = incidentRay - 2 * projectedPoint.intersectedTriangle.normal * (glm::dot(incidentRay, projectedPoint.intersectedTriangle.normal));
	// normalised ray from point to camera (negative of the directionVector)
	float brightness = clamp(glm::dot(-view, reflectedRay), 0, 1);
	// raise brightness value to given power to adjust specular spread
	brightness = pow(brightness, reflectiveness);
	// brightness *= SPECINT; 
	brightness *= thisLight.intensity;
	return brightness;
} 

float specularLighting(RayTriangleIntersection &projectedPoint, glm::vec3 view, glm::vec3 &pointNormal, Light &thisLight, int reflectiveness = SPECPOW) {
	// cast ray from selected light to point
	glm::vec3 incidentRay = glm::normalize(projectedPoint.intersectionPoint - thisLight.light);
	// reflected ray is derived from following formula (Ri - 2 N (Ri . N))
	glm::vec3 reflectedRay = incidentRay - 2 * pointNormal * (glm::dot(incidentRay, pointNormal));
	// normalised ray from point to camera (negative of the directionVector)
	float brightness = clamp(glm::dot(-view, reflectedRay), 0, 1);
	// raise brightness value to given power to adjust specular spread
	brightness = pow(brightness, reflectiveness);
	// if (brightness > 0.9) std::cout << brightness << std::endl;
	// brightness *= SPECINT;
	brightness *= thisLight.intensity;
	return brightness;
} 


float incidenceLighting(RayTriangleIntersection &projectedPoint, Light &thisLight) {
	// cast ray from point to selected light
	glm::vec3 rayToLight = glm::normalize(thisLight.light - projectedPoint.intersectionPoint);
	// raytolight and pointnormal are normalized so dot product returns cos(angle)
	float angle = glm::dot(rayToLight, projectedPoint.intersectedTriangle.normal);
	// clamp min set to above 0 so that surfaces not incident (obtuse angle - facing away) from the light aren't black
	angle *= thisLight.intensity;
	return clamp(angle, 0.1, 1);
}

float incidenceLighting(RayTriangleIntersection &projectedPoint, glm::vec3 &pointNormal, Light &thisLight) {
	// cast ray from point to selected light
	glm::vec3 rayToLight = glm::normalize(thisLight.light - projectedPoint.intersectionPoint);
	// raytolight and pointnormal are normalized so dot product returns cos(angle)
	float angle = glm::dot(rayToLight, pointNormal);
	// clamp min set to above 0 so that surfaces not incident (obtuse angle - facing away) from the light aren't black
	angle *= thisLight.intensity;
	return clamp(angle, 0.1, 1);
}

float proximityLighting(RayTriangleIntersection &projectedPoint, Light &thisLight, float proxintensity = PROXINT) {
	// cast ray from point to selected light
	glm::vec3 rayToLight = thisLight.light - projectedPoint.intersectionPoint;
	// radius of light (distance from light to point)
	float radius = glm::length(rayToLight);

	float brightness = proxintensity / (4 * M_PI * radius * radius);
	// clamp min set to 0 so that VERY far away spots can be nearly black
	brightness *= thisLight.intensity;
	return clamp(brightness, 0, 1);
}

float shadowTrace(RayTriangleIntersection &projectedPoint, Light &thisLight) {
	float shadowVal = 1;
	// cast shadow ray from point to selected light
	glm::vec3 shadowRay = thisLight.light - projectedPoint.intersectionPoint;
	float distanceToLight = glm::length(shadowRay);
	shadowRay = glm::normalize(shadowRay);
	
	if(getClosestIntersection(shadowRay, projectedPoint.intersectionPoint, projectedPoint.triangleIndex).distanceFromCamera < distanceToLight) {
		shadowVal = 1 - 0.5 * thisLight.intensity;
	}
	else {
		if (thisLight.intensity > 0) shadowVal = 1 / (1 - 0.5 * thisLight.intensity);
	}
	
	return clamp(shadowVal,0,1);
	// RayTriangleIntersection shadowIntersection = getClosestIntersection(shadowRay, projectedPoint.intersectionPoint, projectedPoint.triangleIndex);
	// return (shadowIntersection.distanceFromCamera < distanceToLight);
	// return (getClosestIntersection(shadowRay, projectedPoint.intersectionPoint, projectedPoint.triangleIndex).distanceFromCamera < distanceToLight);
}

float gouraudShading(RayTriangleIntersection &projectedPoint, glm::vec3 camToPoint, Colour &extraColour) {
	// use intersection.solution to get barycentric coords of triangle
	float u, v, w;
	u = projectedPoint.solution.y;
	v = projectedPoint.solution.z;
	w = 1 - (u + v);
	

	glm::vec3 pointNormal(0,0,0);
	pointNormal += projectedPoint.intersectedTriangle.vertexNormals[0] * w;
	pointNormal += projectedPoint.intersectedTriangle.vertexNormals[1] * u;
	pointNormal += projectedPoint.intersectedTriangle.vertexNormals[2] * v;
	pointNormal = glm::normalize(pointNormal);

	float proximity, incidence, specular, ambient = 0.1;

	float brightestLight = 0;
	float brightnessCoeff = 1;
	float proxVal = 0;
	float incVal = 0;
	float specVal = 0;


	for (int lt = 0; lt < lightList.size(); lt++) {
		Light thisLight = lightList.at(lt);
		brightestLight = std::max(brightestLight, thisLight.intensity);

		proximity = proximityLighting(projectedPoint, thisLight, 8);
		proxVal = std::max(proxVal, proximity);
		incidence = incidenceLighting(projectedPoint, pointNormal, thisLight);
		incVal = std::max(incVal, incidence);
		specular = specularLighting(projectedPoint, camToPoint, pointNormal, thisLight, 64);
		specVal = std::max(specVal, specular);

		float diffMax = std::max({0.0f, proximity * incidence});
		float specMax = std::max({0.00f, specular});
		float diffWeight = 0.4;
		float specWeight = 1;
		float ambWeight = 0.5 * thisLight.intensity;
		// if (specular > 0.2) std::cout << " sp: " << specular << std::endl;
		// float scale = thisLight.colourIntensity * thisMax;
		// float scale = thisMax * thisLight.intensity;
		// if (scale > 0.5) std::cout <<"scale: " << scale << std::endl;

		extraColour.red += thisLight.colour.red * (diffMax * diffWeight + specMax * specWeight + ambient * ambWeight);
		extraColour.green += thisLight.colour.green * (diffMax * diffWeight + specMax * specWeight + ambient * ambWeight);
		extraColour.blue += thisLight.colour.blue * (diffMax * diffWeight + specMax * specWeight + ambient * ambWeight);

	}

	brightnessCoeff *= brightestLight;
	return clamp(std::max({brightnessCoeff * incVal * proxVal, brightnessCoeff * specVal, brightnessCoeff * ambient}), 0, 1);
	// return clamp(std::max({ambient, specular}), 0, 1);
	// return incidence;
	// return specular;
	// return proximity;
}

void raytraceScene(DrawingWindow &window, bool shadowsOn = false, bool texturesOn = false, bool useProximity = false, bool useIncidence = false, bool useSpecular = false) {
	for (int ys = 0; ys < HEIGHT; ys++) {
		for (int xs = 0; xs < WIDTH; xs++) {
			// get vector from camera to point on plane
			// point on plane defined
			glm::vec3 planePoint(xs -float(WIDTH/2), - ys + HEIGHT/2, -FOCLEN);
			planePoint.x /= (SCALE);
			planePoint.y /= (SCALE);
			// vector from camera position to point
			glm::vec3 directionVector = planePoint - campos;
			directionVector = camrot * directionVector;
			directionVector = glm::normalize(directionVector);

			RayTriangleIntersection intersection = getClosestIntersection(directionVector);
			if (inbounds(xs, ys)) {
				// pixel colour is the colour of the pixel at (xs,ys) and we can alter it depending if it is in shadow
				Colour pixelColour(0,0,0);
				// 	// if (inrange(texturePoints[0].x, texturePoints[0].y, 0, 1)) {

				if (intersection.distanceFromCamera < INFINITY) {
					// if triangle is textured, we need to get the pixel colour from the texture
					if (intersection.intersectedTriangle.textured && texturesOn) {
						// pointer to texturemap in global memory instead of making a new texture every time
						TextureMap *texture = &textureList.at(intersection.intersectedTriangle.textureIndex);
						std::array<TexturePoint,3> *texturePoints = &intersection.intersectedTriangle.texturePoints;

						glm::vec2 e0((*texturePoints)[1].x - (*texturePoints)[0].x, (*texturePoints)[1].y - (*texturePoints)[0].y);
						glm::vec2 e1((*texturePoints)[2].x - (*texturePoints)[0].x, (*texturePoints)[2].y - (*texturePoints)[0].y);
						glm::vec2 texturePixel((*texturePoints)[0].x, (*texturePoints)[0].y);
						texturePixel += (intersection.solution.y * e0 + intersection.solution.z * e1);
						texturePixel.x = clamp(texturePixel.x * (*texture).width, 0, (*texture).width - 1);
						texturePixel.y = clamp(texturePixel.y * (*texture).height, 0, (*texture).height - 1);
						// texturePixel.y *= (*texture).height;
						uint32_t rgbval = (*texture).pixels[(round(texturePixel.y) * (*texture).width) + round(texturePixel.x)];;
						pixelColour = Colour(rgbval);
						// std::cout << " x: " << texturePixel.x <<  " y: " << texturePixel.y << std::endl;
					}
					// else if untextured, set colour to triangle colour
					else pixelColour = intersection.intersectedTriangle.colour;
					
					float brightnessCoeff = 1;
					// min value of lighting 
					float ambientLight = 0.2;
					float specularCoeff = 0;
					// float diffuse = 1;
					// float specular = 0;
					std::vector<Colour> lightColours;
					Colour extraColour(0,0,0);

					// if using gouraud shading
					if (intersection.intersectedTriangle.shaded) {
						auto gouraudResult = gouraudShading(intersection, directionVector, extraColour);
						brightnessCoeff = gouraudResult;
					}
					// else use standard lighting 
					// we need to dim the brightness of colour of the pixels in shadow by using a brightness coefficient
					else {
						// line-of-sight-shadows
						// float shadowVal = 1;
						float brightestLight = 0;
						float proxVal = 0;
						float incVal = 0;
						float specVal = 0;

						for (int lt = 0; lt < lightList.size(); lt++) {
							Light thisLight = lightList.at(lt);
							float thisProx = 0;
							float thisInc = 0;
							float thisSpec = 0;
							float shadowVal = 0;

							if (shadowsOn) {
								shadowVal = shadowTrace(intersection, thisLight);
								brightestLight = std::max(thisLight.intensity, brightestLight);
								brightnessCoeff *= shadowVal;
							}
							// proximity lighting
							if (useProximity) {
								thisProx = proximityLighting(intersection, thisLight);
								proxVal = std::max(thisProx, proxVal);
							}
							// incidence lighting
							if (useIncidence) {
								thisInc = incidenceLighting(intersection, thisLight);
								incVal = std::max(thisInc, incVal);
							}

							if (useSpecular) {
								thisSpec = specularLighting(intersection, directionVector, thisLight);
								specVal = std::max(thisSpec, specVal);
							}

							float thisMax = std::max({ambientLight, thisProx * thisInc, thisSpec});
							float scale = thisLight.colourIntensity * shadowVal * thisMax * thisLight.intensity;
							extraColour.red += thisLight.colour.red * scale;
							extraColour.green += thisLight.colour.green * scale;
							extraColour.blue += thisLight.colour.blue * scale;

						}
						brightnessCoeff *= brightestLight;
						specularCoeff = clamp(specVal * SPECINT, 0, 255);
						brightnessCoeff = std::max({brightnessCoeff * ambientLight, brightnessCoeff * proxVal * incVal, brightnessCoeff * specVal});
						clamp(brightnessCoeff,0,1);
					}

					// white lights 
					pixelColour.red = clamp(pixelColour.red * brightnessCoeff + specularCoeff + extraColour.red, 0, 255);
					pixelColour.green = clamp(pixelColour.green * brightnessCoeff + specularCoeff + extraColour.green, 0, 255);
					pixelColour.blue = clamp(pixelColour.blue * brightnessCoeff + specularCoeff + extraColour.blue, 0, 255);

				}
				window.setPixelColour(xs, ys, ColourToInt(pixelColour));
			}
		}
	}
}

CanvasPoint scalePoint(CanvasPoint p, float scaleFactor) {
	p.x *= scaleFactor;
	p.y *= scaleFactor;
	return p;
}

std::vector<ModelTriangle> scaleModelTriangles(std::vector<ModelTriangle> triangleList, float scaleFactor) {
	for (int i = 0; i < triangleList.size(); i++) {
		triangleList[i].vertices[0] *= scaleFactor;
		triangleList[i].vertices[1] *= scaleFactor;
		triangleList[i].vertices[2] *= scaleFactor;
	}
	return triangleList;
}

void clearScene(DrawingWindow &window) {
	window.clearPixels();
	for (int i = 0; i < WIDTH * HEIGHT; i++) {
		depthBuffer[i] = 0;
	}
}

void RTlookAt(glm::vec3 target) {
	glm::vec3 forward, up, right;
	forward = glm::normalize(campos - target);
	right = glm::normalize(glm::cross(glm::vec3(0,1,0), forward));
	up = glm::normalize(glm::cross(forward, right));

	// std::cout << vecToCP(forward) << vecToCP(right) << vecToCP(up) << std::endl;
	camrot = glm::mat3(right, up, forward);
}

void record(DrawingWindow &window) {
	if (frame == 0 && keyFrame == 0) {

		// posstep = (positions[keyFrame + 1] - positions[keyFrame]) / keyframeTime[keyFrame];
		// rotstep = (rotations[keyFrame + 1] - rotations[keyFrame]) / keyframeTime[keyFrame];
		l1step = (L1Int[keyFrame + 1] - L1Int[keyFrame]) / keyframeTime[keyFrame];
		l2step = (L2Int[keyFrame + 1] - L2Int[keyFrame]) / keyframeTime[keyFrame];
		
		l1colstepr = (L1Col[1].red - L1Col[0].red) / keyframeTime[5];
		l1colstepg = (L1Col[1].green - L1Col[0].green) / keyframeTime[5];
		l1colstepb = (L1Col[1].blue - L1Col[0].blue) / keyframeTime[5];

		l2colstepr = (L2Col[1].red - L2Col[0].red) / keyframeTime[5];
		l2colstepg = (L2Col[1].green - L2Col[0].green) / keyframeTime[5];
		l2colstepb = (L2Col[1].blue - L2Col[0].blue) / keyframeTime[5];

	}

	// else if (keyFrame > keyframeTime.size()) {
	// 	return;
	// }

	else if (frame >= keyframeTime[keyFrame]) {
		keyFrame++;
		frame = 0;
		// posstep = (positions[keyFrame + 1] - positions[keyFrame]) / keyframeTime[keyFrame];
		// rotstep = (rotations[keyFrame + 1] - rotations[keyFrame]) / keyframeTime[keyFrame];
		l1step = (L1Int[keyFrame + 1] - L1Int[keyFrame]) / keyframeTime[keyFrame];
		l2step = (L2Int[keyFrame + 1] - L2Int[keyFrame]) / keyframeTime[keyFrame];
	}

	if (keyFrame == 5) {
		lightList[0].intensity = 0.8;
		lightList[1].intensity = 0.8;
		lightList[0].colourIntensity = 2;
		lightList[1].colourIntensity = 1;

		lightList[0].colour.red = L1Col[0].red + l1colstepr * frame;
		lightList[0].colour.green = L1Col[0].green + l1colstepg * frame;
		lightList[0].colour.blue = L1Col[0].blue + l1colstepb * frame;

		lightList[1].colour.red = L2Col[0].red + l2colstepr * frame;
		lightList[1].colour.green = L2Col[0].green + l2colstepg * frame;
		lightList[1].colour.blue = L2Col[0].blue + l2colstepb * frame;

	}
	else {
		lightList[0].intensity = L1Int[keyFrame] + (l1step * frame);
		lightList[1].intensity = L2Int[keyFrame] + (l2step * frame);
	}
	// campos = positions[keyFrame] + (posstep * frame);
	// camrot = rotations[keyFrame] + (rotstep * frame);

	std::string filename = "./lighting/" + std::to_string(keyFrame) + '_' + std::to_string(frame);
	window.saveBMP(filename);
	std::cout << filename << std::endl;
	frame++;
}

void draw(DrawingWindow &window) {

	// if rendering scene as a whole (raytracing)
	if (renderStyle == 3) {
		raytraceScene(window);
	}
	// raytracing with all on
	else if (renderStyle == 4) {

		// // snippet to make ball reflective
		// for (int i = 32; i < 144; i++) {
		// 	triangleList[i].reflectiveness = 1;
		// }
		raytraceScene(window, true, true, true, true, true);
	}
	// raytracing with los shadows, no texture, proximity lighting
	else if (renderStyle == 5) {
		raytraceScene(window, true, true, true, true, true);
		if (keyFrame < keyframeTime.size())	record(window);
	}
	// raytracing with no los shadows, no texture, no proximity lighting, incidence lighting
	else if (renderStyle == 6) {
		raytraceScene(window, false, false, false, false, true);
	}
	else if (renderStyle == 7) {
		raytraceScene(window, true, true, true, true, true);
	}
	// else looping through triangles (rasterising)
	else {
	// for (int i = 0; i < 12; i++) { // no boxes
	// for (int i = 12; i < 22; i++) { // red box
	// for (int i = 22; i < 32; i++) { // blue box
	// for (int i = 12; i < triangleList.size(); i++) { // both boxes

		// code to show lights when wireframed
		triangleList.resize(modelCount);
		for (int i = 0; i < lightList.size(); i++) {
			ModelTriangle lightTriangle;	
			float lx, ly, lz, delta = 0.03, lightint, colint;
			lx = lightList[i].light.x;
			ly = lightList[i].light.y;
			lz = lightList[i].light.z;
			lightint = lightList[i].intensity;
			colint = lightList[i].colourIntensity;
			glm::vec3 a(lx - delta, ly, lz), b(lx, ly + delta, lz), c(lx + delta, ly, lz);

			lightTriangle.vertices[0] = a;
			lightTriangle.vertices[1] = b;
			lightTriangle.vertices[2] = c;
			lightTriangle.colour = lightList[i].colour;
			lightTriangle.colour.red *= lightint;
			lightTriangle.colour.green *= lightint;
			lightTriangle.colour.blue *= lightint;
			
			triangleList.push_back(lightTriangle);
		}

		for (int i = 0; i < triangleList.size(); i++) { //whole scene
			ModelTriangle t = triangleList[i];
			// CanvasTriangle t;
			Colour colour(triangleList[i].colour);

			for (int j = 0; j < 3; j++) {
				CanvasPoint p = getCanvasIntersectionPoint(triangleList[i].vertices[j]);
				t.vertices[j].x = p.x;
				t.vertices[j].y = p.y;
				t.vertices[j].z = p.depth;
			}

			// switch case for changing render style 
			// 1 - wireframes
			// 2 - rasterising 
			// 3 - raytracing
			// std::string colourname = colour.name;
			switch (renderStyle) {
				case 1:
					drawStrokedTriangle(t, window, colour);
					break;
				
				case 2: 
					/// NEED TO IMPLEMENT TEXTURE DETECTION + GET TEXTURES DRAWING PROPERLY
					if (t.texturePoints[0].x > 0) {
						textureTriangle(t, window);
					}
					else {
						rasterizeTriangle(modelToCanvas(t), window, colour);
					}
					break;

				default:
					break;
			}
		}
		triangleList.resize(modelCount);
	}
}

void readMaterialFile(std::unordered_map<std::string, Colour> &materials, std::string filename) {
	std::ifstream file(filename, std::ifstream::in);	
	std::string readLine;
	std::string name;
	Colour colour;
	// int textureCount = 0;
	textureList.clear();

	if (file.is_open()) {
		while(std::getline(file, readLine)) {
			std::vector<std::string> lineSegments = split(readLine,' ');
			if (lineSegments[0] == "newmtl") {
				name = lineSegments[1];
			}
			else if (lineSegments[0] == "Kd") {
				int r, g, b;
				r = std::stof(lineSegments[1]) * 255;
				g = std::stof(lineSegments[2]) * 255;
				b = std::stof(lineSegments[3]) * 255;
				colour = Colour(name,r,g,b);
				materials.emplace(name,colour);
			}
			else if (lineSegments[0] == "map_Kd") {
				// add the texture to the textureList
				std::string texturePath = MODELS;
				texturePath.append("/" + lineSegments[1]);
				textureList.emplace_back(TextureMap(texturePath));
				
				// textureCount++;
			}
		}
	}
} 

void readOBJFile(std::string filename, DrawingWindow &window, std::unordered_map<std::string, Colour> materials, glm::vec3 transpose = glm::vec3(0, 0, 0), float scale = 1, int currentTextureIndex = 0) {
	std::ifstream file(filename, std::ifstream::in);	
	std::string readLine;
	std::vector<glm::vec3> vertexList;
	std::vector<glm::vec3> vertexNormalList;
	std::vector<TexturePoint> texturePointList;
	Colour currentMaterial(255,0,0);
	bool currentMaterialIsTexture = false;

	if (file.is_open()) {
		while(std::getline(file, readLine)) {
			std::vector<std::string> lineSegments = split(readLine, ' ');

			// if (lineSegments[0] == "mtllib") {
			// 	std::string mtllib = filename;
			// 	int folderindex = mtllib.find_last_of('/');
			// 	mtllib.erase(folderindex, mtllib.size());
			// 	mtllib.append("/" + lineSegments[1]);
			// 	readMaterialFile(materials, mtllib);
			// }
			if (lineSegments[0] == "v") {
				float x, y, z;
				x = std::stof(lineSegments[1]);
				y = std::stof(lineSegments[2]);
				z = std::stof(lineSegments[3]);

				vertexList.push_back(glm::vec3(x,y,z)*scale + transpose);
			}
			else if (lineSegments[0] == "vt") {
				float x, y;
				x = std::stof(lineSegments[1]);
				y = std::stof(lineSegments[2]);
				texturePointList.push_back(TexturePoint(x,y));
				currentMaterialIsTexture = true;
			}
			else if (lineSegments[0] == "vn") {
				float x, y, z;
				x = std::stof(lineSegments[1]);
				y = std::stof(lineSegments[2]);
				z = std::stof(lineSegments[3]);

				vertexNormalList.push_back(glm::vec3(x,y,z));
			}

			else if (lineSegments[0] == "f") {
				int a, b, c, ta, tb, tc;
				a = std::stoi(lineSegments[1]) - 1;
				b = std::stoi(lineSegments[2]) - 1;
				c = std::stoi(lineSegments[3]) - 1;

				ModelTriangle face(vertexList[a], vertexList[b], vertexList[c], currentMaterial);

				// checking face for texturepoints
				if (!split(lineSegments[1], '/').at(1).empty()){
					ta = std::stoi(split(lineSegments[1], '/').at(1)) - 1;
					tb = std::stoi(split(lineSegments[2], '/').at(1)) - 1;
					tc = std::stoi(split(lineSegments[3], '/').at(1)) - 1;
					face.texturePoints[0] = texturePointList[ta];
					face.texturePoints[1] = texturePointList[tb];
					face.texturePoints[2] = texturePointList[tc];
					face.textured = true;
					face.textureIndex = currentTextureIndex;
				}
				else face.textured = false;

				// checking face for vertexNormals
				if (split(lineSegments[1], '/').size() == 3) {
					for (int i = 0; i < 3; i++) {
						float vnIndex = stoi(split(lineSegments[i+1], '/').at(2));
						face.vertexNormals[i] = vertexNormalList[vnIndex -1];
					}
					face.shaded = true;
				}
				else face.shaded = false;
				
				// calculate normal of triangle
				// edge AB
				glm::vec3 e0 = vertexList[b] - vertexList[a];
				// edge AC
				glm::vec3 e1 = vertexList[c] - vertexList[a];
				// if vertices are wound clockwise then 
				if (a < b) {
					face.normal = glm::normalize(glm::cross(e0, e1));
				}
				// else vertices are wound anticlockwise 
				else {
					face.normal = glm::normalize(glm::cross(e1, e0));
				}
				triangleList.push_back(face);
			}
			else if (lineSegments[0] == "usemtl") {
				if (currentMaterialIsTexture) currentTextureIndex++;
				currentMaterialIsTexture = false;
				currentMaterial = materials.at(lineSegments[1]);
			}
		}

	}
	else std::cout << "file not open" << std::endl;

	// triangleList.push_back(ModelTriangle(lightList[0]+glm::vec3(0.1,0,0),lightList[0],lightList[0]+glm::vec3(0,-0.1,0),Colour(255,255,255)));

	// triangleList = scaleModelTriangles(triangleList, 2);
}

void loadObjs(DrawingWindow &window) {
	std::unordered_map<std::string, Colour> materials;
	readMaterialFile(materials, MATS);

	triangleList.clear();
	readOBJFile(BOX, window, materials);
	// top of box glm::vec3(-1,0.1,-1.6)
	readOBJFile(BALL, window, materials, glm::vec3(-1.5,-3.3,0.5));
	readOBJFile(LOGO, window, materials, glm::vec3(1,0,1), 0.1, 1);
	triangleList[26].reflectiveness = 1;
	triangleList[31].reflectiveness = 1;
	modelCount = triangleList.size();
	triangleList = scaleModelTriangles(triangleList, 0.17);
}

void lookAt(CanvasPoint p) {
	glm::vec3 forward, up, right, target(p.x, p.y, p.depth);
	forward = glm::normalize(campos - target);
	right = glm::normalize(glm::cross(glm::vec3(0,1,0), forward));
	up = glm::normalize(glm::cross(forward, right));

	// std::cout << vecToCP(forward) << vecToCP(right) << vecToCP(up) << std::endl;
	camrot = glm::mat3(right, up, forward);
}

void rotateScene(int dir = 0, float angle = M_PI/16) {
	glm::mat3 matrix;
	switch (dir) {
	case 1:
		// rotate left (about y) key: J
		matrix = glm::mat3(
			cos(angle), 0, -sin(angle),
			         0, 1, 0,
			sin(angle), 0, cos(angle)
		);
		campos = matrix * campos;
		// camrot = matrix * camrot;
		break;
		
	
	case 2:
		// rotate right (about y) key: L
		matrix = glm::mat3(
			cos(-angle), 0, -sin(-angle),
			          0, 1, 0,
			sin(-angle), 0, cos(-angle)
		);
		campos = matrix * campos;
		// camrot = matrix * camrot;
		break;
	
	case 3:
		// rotate up (about x) key I
		matrix = glm::mat3(
			1,           0, 0,
			0,  cos(angle), sin(angle),
			0, -sin(angle), cos(angle)
		);
		campos = matrix * campos;
		// camrot = matrix * camrot;
		break;
	
	case 4:
		// rotate down (about x) key K
		matrix = glm::mat3(
			1,            0, 0,
			0,  cos(-angle), sin(-angle),
			0, -sin(-angle), cos(-angle)
		);		
		campos = matrix * campos;
		// camrot = matrix * camrot;
		break;
	
	default:
		break;
	}
	lookAt(ORIGIN);
}

void orbitCamera(DrawingWindow &window, int f) {
	int delta = 64;
	rotateScene(1, 2 * M_PI / delta);
	lookAt(ORIGIN);
	
	if (f >= delta) renderStyle = 2;

	std::string filename = "./orbit/" + std::to_string(f);
	window.saveBMP(filename);
	std::cout << filename << std::endl;

	if (f > 2 * delta) OrbiterToggle = false;
}

void nudgeCamera(int dir = 0, float mag = 0.1) {
	switch (dir)
	{
	case 1:
		// move camera left
		campos.x -= mag;
		break;
	
	case 2:
		// move camera right
		campos.x += mag;
		break;

	case 3:
		// move camera up
		campos.y -= mag;
		break;

	case 4:
		// move camera down
		campos.y += mag;
		break;

	case 5: 
		// forward
		campos.z -= 3*mag;
		break;

	case 6:
		// backwards
		campos.z += 3*mag;
		break;

	default:
		break;
	}
}

void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		// Move camera in XY with arrow keys
		if (event.key.keysym.sym == SDLK_LEFT) {
			nudgeCamera(1);
		}
		else if (event.key.keysym.sym == SDLK_RIGHT) {
			nudgeCamera(2);
		}
		else if (event.key.keysym.sym == SDLK_UP) {
			nudgeCamera(4);
		}
		else if (event.key.keysym.sym == SDLK_DOWN) {
			nudgeCamera(3);
		}
		else if (event.key.keysym.sym == SDLK_w) {
			nudgeCamera(5);
		}
		else if (event.key.keysym.sym == SDLK_s) {
			nudgeCamera(6);
		}

		// Rotate camera in XY with IJKL
		else if (event.key.keysym.sym == SDLK_j ) {
			rotateScene(1);
		}
		else if (event.key.keysym.sym == SDLK_l) {
			rotateScene(2);
		}
		else if (event.key.keysym.sym == SDLK_i ) {
			rotateScene(3);
		}
		else if (event.key.keysym.sym == SDLK_k) {
			rotateScene(4);
		}
		// Camera orbiter toggle Key = P
		else if (event.key.keysym.sym == SDLK_p) {
			OrbiterToggle = !OrbiterToggle;
		}
		else if (event.key.keysym.sym == SDLK_SEMICOLON) {
			lookAt(CanvasPoint(0,0,0));
		}

		// change render style with numbers
		// 1 = wireframes
		else if (event.key.keysym.sym == SDLK_1) {
			renderStyle = 1;
		}
		// 2 = rasterise
		else if (event.key.keysym.sym == SDLK_2) {
			renderStyle = 2;
		}
		else if (event.key.keysym.sym == SDLK_3) {
			renderStyle = 3;
		}
		else if (event.key.keysym.sym == SDLK_4) {
			renderStyle = 4;
		}
		else if (event.key.keysym.sym == SDLK_5) {
			renderStyle = 5;
		}
		else if (event.key.keysym.sym == SDLK_6) {
			renderStyle = 6;
		}
		else if (event.key.keysym.sym == SDLK_7) {
			renderStyle = 7;
		}
		else if (event.key.keysym.sym == SDLK_8) {
			renderStyle = 8;
		}
		else if (event.key.keysym.sym == SDLK_9) {
			renderStyle = 9;
		}
		else if (event.key.keysym.sym == SDLK_0) {
			renderStyle = 0;
		}

		else if (event.key.keysym.sym == SDLK_KP_6) lightList[0].light.x += 0.1;
		else if (event.key.keysym.sym == SDLK_KP_4) lightList[0].light.x -= 0.1;
		else if (event.key.keysym.sym == SDLK_KP_8) lightList[0].light.y += 0.1;
		else if (event.key.keysym.sym == SDLK_KP_5) lightList[0].light.y -= 0.1;
		else if (event.key.keysym.sym == SDLK_KP_7) lightList[0].light.z += 0.1;
		else if (event.key.keysym.sym == SDLK_KP_1) lightList[0].light.z -= 0.1;
		else if (event.key.keysym.sym == SDLK_KP_9) lightList[1].intensity = clamp(lightList[1].intensity + 0.1, 0, 1);
		else if (event.key.keysym.sym == SDLK_KP_3) lightList[1].intensity = clamp(lightList[1].intensity - 0.1, 0, 1);
		else if (event.key.keysym.sym == SDLK_KP_PLUS) lightList[0].intensity = clamp(lightList[0].intensity + 0.1, 0, 1);
		else if (event.key.keysym.sym == SDLK_KP_ENTER) lightList[0].intensity = clamp(lightList[0].intensity - 0.1, 0, 1);

		else if (event.key.keysym.sym == SDLK_z) {
			std::cout << "campos: " << "(" << campos.x << ", " << campos.y << ", " << campos.z << ")" << std::endl;
			std::cout << "camrot: ("; 
			for (int i = 0; i < 3; i++) {
				glm::vec3 v = glm::column(camrot, i);
				std::cout << v.x << ", " << v.y << ", " << v.z << ", ";
			}
			std::cout << ")" << std::endl; 
		}

		else if (event.key.keysym.sym == SDLK_x) {
			keyFrame++;
		}


		// else if (event.key.keysym.sym == SDLK_m) readMaterialFile(std::unordered_map<std::string, Colour> materials, MATS);
		else if (event.key.keysym.sym == SDLK_o) {
			loadObjs(window);
			draw(window);
		} 
		// else if (event.key.keysym.sym == SDLK_c) clearScene(window);
		else if (event.key.keysym.sym == SDLK_r) {
			clearScene(window);
			triangleList.clear();
			campos = DEFPOS;
			camrot = DEFROT;
			lightList = DEFLIGHT;
		}

	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	int f = 0;

	// std::cout << BRICK << BOX << MATS << std::endl;
	loadObjs(window);
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		if (OrbiterToggle) {
			orbitCamera(window, f);
			f++;
		}	
		clearScene(window);
		draw(window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}