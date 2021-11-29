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
#include <fstream>
#include <string>
#include <iostream>
#include <unordered_map>
#include <cmath>

#define WIDTH 300
#define HEIGHT 300
// #define BRICK "../textures/texture.ppm"
#define MODELS "/home/jakehughes/Desktop/CGCW/models"
#define BOX "/home/jakehughes/Desktop/CGCW/models/textured-cornell-box.obj"
// #define BOX "/home/jakehughes/Desktop/CG2021/Weekly Workbooks/04 Wireframes and Rasterising/Wireframes/models/boxes.obj"
#define MATS "/home/jakehughes/Desktop/CGCW/models/textured-cornell-box.mtl"
#define FOCLEN 6 // default focal length
#define SCALE 50
#define DEFPOS glm::vec3(0.0, 0.0, 10.0)
#define DEFROT glm::mat3(1,0,0,0,1,0,0,0,1)
#define ORIGIN CanvasPoint(0,0,0)

// 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

bool OrbiterToggle;
glm::vec3 campos(0.0, 0.0, 10.0); 
float depthBuffer[WIDTH * HEIGHT];
std::vector<ModelTriangle> triangleList;
std::vector<std::string> textureList;
glm::mat3 camrot(DEFROT);
int renderStyle = 0;

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

uint32_t getTexturePixel(TexturePoint p, TextureMap &map) {
	return map.pixels[(round(p.y) * map.width) + round(p.x)];
}

void textureTriangle(ModelTriangle t, DrawingWindow &window, int textureIndex = 0) {
	// std::cout << textureList.at(textureIndex) << std::endl;
	TextureMap texture(textureList.at(textureIndex));
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
		topT.texturePoints = {intercept.texturePoint, t.texturePoints[0], t.texturePoints[1]};
		topT = sortTopTriangle(topT);

		botT = ModelTriangle(cpToVec(intercept), t.vertices[1], t.vertices[2], t.colour);
		botT.texturePoints = {intercept.texturePoint, t.texturePoints[1], t.texturePoints[2]};
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

CanvasPoint scalePoint(CanvasPoint p, float scaleFactor) {
	p.x *= scaleFactor;
	p.y *= scaleFactor;
	return p;
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

void draw(DrawingWindow &window) {
	// for (int i = 0; i < 12; i++) { // no boxes
	// for (int i = 12; i < 22; i++) { // red box
	// for (int i = 22; i < 32; i++) { // blue box
	// for (int i = 12; i < triangleList.size(); i++) { // both boxes
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
}

void readMaterialFile(std::unordered_map<std::string, Colour> &materials, std::string filename) {
	std::ifstream file(filename, std::ifstream::in);	
	std::string readLine;
	std::string name;
	Colour colour;
	int textureCount = 0;
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
				textureList.emplace_back(texturePath);
				
				textureCount++;
			}
		}
	}
} 

void readOBJFile(std::string filename, DrawingWindow &window) {
	std::ifstream file(filename, std::ifstream::in);	
	std::string readLine;
	std::vector<glm::vec3> vertexList;
	std::vector<TexturePoint> texturePointList;
	std::unordered_map<std::string, Colour> materials;
	Colour currentMaterial;

	triangleList.clear();

	if (file.is_open()) {
		while(std::getline(file, readLine)) {
			std::vector<std::string> lineSegments = split(readLine, ' ');

			if (lineSegments[0] == "mtllib") {
				std::string mtllib = filename;
				int folderindex = mtllib.find_last_of('/');
				mtllib.erase(folderindex, mtllib.size());
				mtllib.append("/" + lineSegments[1]);
				readMaterialFile(materials, mtllib);
			}
			if (lineSegments[0] == "v") {
				float x, y, z;
				x = std::stof(lineSegments[1]);
				y = std::stof(lineSegments[2]);
				z = std::stof(lineSegments[3]);

				vertexList.push_back(glm::vec3(x,y,z));
			}
			else if (lineSegments[0] == "vt") {
				float x, y;
				x = std::stof(lineSegments[1]);
				y = std::stof(lineSegments[2]);
				texturePointList.push_back(TexturePoint(x,y));
			}

			else if (lineSegments[0] == "f") {
				int a, b, c, ta, tb, tc;
				a = std::stoi(lineSegments[1]) - 1;
				b = std::stoi(lineSegments[2]) - 1;
				c = std::stoi(lineSegments[3]) - 1;

				ModelTriangle face(vertexList[a], vertexList[b], vertexList[c], currentMaterial);

				if (!split(lineSegments[1], '/').at(1).empty()){
					ta = std::stoi(split(lineSegments[1], '/').at(1)) - 1;
					tb = std::stoi(split(lineSegments[2], '/').at(1)) - 1;
					tc = std::stoi(split(lineSegments[3], '/').at(1)) - 1;
					face.texturePoints[0] = texturePointList[ta];
					face.texturePoints[1] = texturePointList[tb];
					face.texturePoints[2] = texturePointList[tc];
				}
				else {
					TexturePoint tp(-1,-1);
					face.texturePoints[0] = tp;
					face.texturePoints[1] = tp;
					face.texturePoints[2] = tp;

				}	

				triangleList.push_back(face);
			}
			else if (lineSegments[0] == "usemtl") {
				currentMaterial = materials.at(lineSegments[1]);
			}
		}

	}
	else std::cout << "file not open" << std::endl;

	triangleList = scaleModelTriangles(triangleList, 0.17);
}

void lookAt(CanvasPoint p) {
	glm::vec3 forward, up, right, target;
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

void rotateCamera(int dir = 0, float angle = M_PI/16) {
	glm::mat3 matrix;
	switch (dir) {
	case 1:
		// rotate left (about y) key: J
		matrix = glm::mat3(
			cos(angle), 0, -sin(angle),
			         0, 1, 0,
			sin(angle), 0, cos(angle)
		);
		// campos = matrix * campos;
		camrot = matrix * camrot;
		break;
		
	
	case 2:
		// rotate right (about y) key: L
		matrix = glm::mat3(
			cos(-angle), 0, -sin(-angle),
			          0, 1, 0,
			sin(-angle), 0, cos(-angle)
		);
		// campos = matrix * campos;
		camrot = matrix * camrot;
		break;
	
	case 3:
		// rotate up (about x) key I
		matrix = glm::mat3(
			1,           0, 0,
			0,  cos(angle), sin(angle),
			0, -sin(angle), cos(angle)
		);
		// campos = matrix * campos;
		camrot = matrix * camrot;
		break;
	
	case 4:
		// rotate down (about x) key K
		matrix = glm::mat3(
			1,            0, 0,
			0,  cos(-angle), sin(-angle),
			0, -sin(-angle), cos(-angle)
		);		
		// campos = matrix * campos;
		camrot = matrix * camrot;
		break;
	
	default:
		break;
	}
}

void orbitCamera() {
	rotateScene(1, 0.01);
	lookAt(ORIGIN);
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
		else if (event.key.keysym.sym == SDLK_0) {
			renderStyle = 0;
		}

		// else if (event.key.keysym.sym == SDLK_m) readMaterialFile(std::unordered_map<std::string, Colour> materials, MATS);
		else if (event.key.keysym.sym == SDLK_o) {
			readOBJFile(BOX, window);
			draw(window);
		} 
		// else if (event.key.keysym.sym == SDLK_c) clearScene(window);
		else if (event.key.keysym.sym == SDLK_r) {
			clearScene(window);
			triangleList.clear();
			campos = DEFPOS;
			camrot = DEFROT;
		}

	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;

	// std::cout << BRICK << BOX << MATS << std::endl;
	readOBJFile(BOX, window);
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		if (OrbiterToggle) {
			orbitCamera();
		}	
		clearScene(window);
		draw(window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}