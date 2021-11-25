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

#define WIDTH 800
#define HEIGHT 800
#define BRICK "../textures/texture.ppm"
#define BOX "/home/jakehughes/Desktop/CG2021/Weekly Workbooks/04 Wireframes and Rasterising/Wireframes/models/cornell-box.obj"
// #define BOX "/home/jakehughes/Desktop/CG2021/Weekly Workbooks/04 Wireframes and Rasterising/Wireframes/models/boxes.obj"
#define MATS "/home/jakehughes/Desktop/CG2021/Weekly Workbooks/04 Wireframes and Rasterising/Wireframes/models/cornell-box.mtl"
#define FOCLEN 6 // default focal length
#define SCALE 150
#define DEFPOS glm::vec3(0.0, 0.0, 10.0)
#define DEFROT glm::mat3(1,0,0,0,1,0,0,0,1)

// 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

glm::vec3 campos(0.0, 0.0, 10.0); 
float depthBuffer[WIDTH * HEIGHT];
std::vector<ModelTriangle> triangleList;
glm::mat3 camrot(1,0,0,
				 0,1,0,
				 0,0,1);


CanvasPoint vecToCP(glm::vec3 vec) {
	CanvasPoint p(vec.x, vec.y, vec.z);
	return p;
}

glm::vec3 cpToVec(CanvasPoint p) {
	glm::vec3 vec(p.x, p.y, p.depth);
	return vec;
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
		if (interpolatedList.back().x > WIDTH || interpolatedList.back().x < 0 || interpolatedList.back().y < 0 || interpolatedList.back().y > HEIGHT) {
			// std::cout << "x: " << interpolatedList.back().x << " y: " << interpolatedList.back().y << " z: " << interpolatedList.back().z << std::endl;
			// i = numberOfValues;
			break;
		}
		interpolatedList.push_back(interpolatedList.back() + step);
	}
	return interpolatedList;
}

void drawLine(CanvasPoint from, CanvasPoint to, DrawingWindow &window, Colour c = Colour(230,10,230)) {
	int xDiff = round(abs(to.x - from.x)); 
	int yDiff = round(abs(to.y - from.y));
	// int zDiff = round(abs(to.depth - from.depth));
	int numberOfValues = std::max(xDiff, yDiff);
	// numberOfValues = std::max(numberOfValues, zDiff);
	numberOfValues *= 2;
	uint32_t colour = ColourToInt(c);

	glm::vec3 vecfrom = cpToVec(from);
	glm::vec3 vecto = cpToVec(to);

	// if both points are on canvas then draw normally
	if (inbounds(from) && inbounds(to)) {
		std::vector<glm::vec3> Points = interpolateThreeElementValues(vecfrom, vecto, numberOfValues);

		// for (int i = 0; i < numberOfValues - 1; i++) {
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
}

void drawStrokedTriangle(CanvasTriangle t, DrawingWindow &window, Colour c = Colour(230,10,230)) {
	drawLine(t.v0(), t.v1(), window, c);
	drawLine(t.v0(), t.v2(), window, c);
	drawLine(t.v1(), t.v2(), window, c);
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

	// //test vertex sort function
	// std::cout << t << std::endl;
	// assert(t.v0().y<=t.v1().y);
	// assert(t.v0().y<=t.v2().y);
	// assert(t.v1().y<=t.v2().y);

	return t;
}

CanvasTriangle sortTopTriangle(CanvasTriangle t) {
	t = sortVerticesVertically(t);
	if(t.v1().x > t.v2().x) std::swap(t.v1(),t.v2());
	return t;
}

CanvasTriangle sortBottomTriangle(CanvasTriangle t) {
	t = sortVerticesVertically(t);
	std::swap(t.v2(), t.v0());
	if(t.v1().x > t.v2().x) std::swap(t.v1(),t.v2());
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

CanvasPoint findIntercept(CanvasTriangle t) {
	// float dx, dy, dix, diy;
	// dx = end.x - start.x;
	// dy = end.y - start.y;
	// diy = start.y - yValue;
	// dix = diy ;

	// if (abs(dy) <= 1) return start.x;
	// else if (abs(dx) <= 1) return start.x + (dx/2);
	// else return start.x + dix;
	
	float verticalHeight = (abs(t.v2().y - t.v0().y));
	float interceptHeight = (abs(t.v1().y - t.v0().y));
	float ratio = interceptHeight / verticalHeight;

	CanvasPoint intercept(lineRatio(t.v0(), t.v2(), ratio));
	return intercept;
}

/* float findXIntercept(float yValue, TexturePoint start, TexturePoint end) {
	float intercept = (yValue * (start.x - end.x) +  ((end.x * start.y) - (start.x * end.y))) / (start.y - end.y);
	return intercept;
} */

/* float findYIntercept(float xValue, CanvasPoint start, CanvasPoint end) {
	// float dx, dy, dix, diy;
	// dx = end.x - start.x;
	// dy = end.y - start.y;
	// dix = start.x - xValue;
	// diy = dix * dy / dx;
	
	// if (abs(dy) <= 1) return start.y;
	// else if (abs(dx) <= 1) return start.y + (dy / 2);
	// else return start.y + diy;
	
	

	float intercept;
	if (abs(start.x - end.x) < 1) {
		intercept = start.y;
	}
	else {
		intercept = (xValue * (start.y - end.y) +  ((start.x * end.y) - (end.x * start.y))) / (start.x - end.x);
	}
	return intercept;
} */

/* float findYIntercept(float xValue, TexturePoint start, TexturePoint end) {
	float intercept = (xValue * (start.y - end.y) +  ((start.x * end.y) - (end.x * start.y))) / (start.x - end.x);
	return intercept;
} */

void rasterizeTriangle(CanvasTriangle t, DrawingWindow &window, Colour c) {

	t = sortVerticesVertically(t); 
	drawStrokedTriangle(t, window, c);

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


	// float minX = std::min(topT.v1().x, topT.v2().x);
	// minX = floor(round(std::min(topT.v0().x, minX)));
	// float maxX = std::max(topT.v1().x, topT.v2().x);
	// maxX = ceil(round(std::max(topT.v0().x, maxX)));
		
	//rasterise top triangle, if there is one
	if (!flatTop) {

		int height = abs(round(topT.v0().y) - round(topT.v1().y)) * 2;
		int base = abs(round(topT.v1().x) - round(topT.v2().x));
		int left = abs(round(topT.v1().x) - round(topT.v0().x));
		int right = abs(round(topT.v0().x - round(topT.v2().x)));

		int width = std::max(left, right);
		width = std::max(width, base);
	

		// std::vector<glm::lowp_vec3>  ApexLeftY = interpolateLowPVec(cpToVec(topT.v0()), cpToVec(topT.v1()), height);
		// std::vector<glm::lowp_vec3> ApexRightY = interpolateLowPVec(cpToVec(topT.v0()), cpToVec(topT.v2()), height);

		// std::vector<glm::lowp_vec3>  LeftApexX = interpolateLowPVec(cpToVec(topT.v1()), cpToVec(topT.v0()), left);
		// std::vector<glm::lowp_vec3> ApexRightX = interpolateLowPVec(cpToVec(topT.v0()), cpToVec(topT.v2()), right);
		// std::vector<glm::lowp_vec3> LeftRight = interpolateLowPVec(cpToVec(topT.v1()), cpToVec(topT.v2()), base);

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

		// std::cout << triangleList[i] << std::endl;
	}
	return triangleList;
}

// std::unordered_map<std::string, Colour> 
void readMaterialFile(std::unordered_map<std::string, Colour> &materials, std::string filename) {
	// std::unordered_map<std::string, Colour> materials;
	std::ifstream file(filename, std::ifstream::in);	
	std::string readLine;
	std::string name;
	Colour colour;


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
		}
	}
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
		CanvasTriangle t;
		for (int j = 0; j < 3; j++) {
			CanvasPoint p = getCanvasIntersectionPoint(triangleList[i].vertices[j]);
			t.vertices[j].x = p.x;
			t.vertices[j].y = p.y;
			t.vertices[j].depth = p.depth;
		}
		// drawStrokedTriangle(t, window, triangleList[i].colour);
		rasterizeTriangle(t, window, triangleList[i].colour);
		// std::cout <<"Done " << i + 1<< " out of " << triangleList.size() << std::endl;
	}
		// std::cout <<"Done All" << std::endl;
}

void readOBJFile(std::string filename, DrawingWindow &window) {
	std::ifstream file(filename, std::ifstream::in);	
	std::string readLine;
	std::vector<glm::vec3> vertexList;
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
			else if (lineSegments[0] == "f") {
				int a, b, c;
				a = std::stoi(lineSegments[1]) - 1;
				b = std::stoi(lineSegments[2]) - 1;
				c = std::stoi(lineSegments[3]) - 1;

				triangleList.push_back(ModelTriangle(vertexList[a], vertexList[b], vertexList[c], currentMaterial));
			}
			else if (lineSegments[0] == "usemtl") {
				// std::cout << '.' << lineSegments[1] << '.' << std::endl;
				currentMaterial = materials.at(lineSegments[1]);
			}
		}

	}
	else std::cout << "file not open" << std::endl;

	triangleList = scaleModelTriangles(triangleList, 0.17);
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

void nudgeCamera(int dir = 0, float mag = 0.1) {
	switch (dir)
	{
	case 1:
		// left
		campos.x += mag;
		break;
	
	case 2:
		// right
		campos.x -= mag;
		break;

	case 3:
		// up
		campos.y += mag;
		break;

	case 4:
		// down
		campos.y -= mag;
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
			rotateCamera(1);
		}
		else if (event.key.keysym.sym == SDLK_l) {
			rotateCamera(2);
		}
		else if (event.key.keysym.sym == SDLK_i ) {
			rotateCamera(3);
		}
		else if (event.key.keysym.sym == SDLK_k) {
			rotateCamera(4);
		}

		// else if (event.key.keysym.sym == SDLK_m) readMaterialFile(std::unordered_map<std::string, Colour> materials, MATS);
		else if (event.key.keysym.sym == SDLK_o) readOBJFile(BOX, window);
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

	std::cout << BRICK << BOX << MATS << std::endl;
	readOBJFile(BOX, window);
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		clearScene(window);
		draw(window);
		window.renderFrame();
	}
}