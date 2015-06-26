// BasicGeometry.cpp
// https://www.topcoder.com/community/data-science/data-science-tutorials/geometry-concepts-basic-concepts/
// https://www.topcoder.com/community/data-science/data-science-tutorials/geometry-concepts-line-intersection-and-its-applications/
// https://www.topcoder.com/community/data-science/data-science-tutorials/geometry-concepts-using-geometry-in-topcoder-problems/


#include "stdafx.h"
#include <utility>
#include <math.h>
#include <vector>
#include <string>

using namespace std;

class Point {
public:
	int x_;
	int y_;
	Point() : x_(0), y_(0) {  }
	Point(int x, int y) : x_(x), y_(y) {}
	Point& operator+(const Point& v) const {
		return Point(x_ + v.x_, y_ + v.y_);
	}
	Point& operator-(const Point&  v) const {
		return Point(x_ + v.x_, y_ + v.y_);
	}
	double magnitude(void) const {
		return sqrt(static_cast<double>(x_ * x_ + y_ * y_));
	}
	int dot_product(const Point& v) const {
		return (x_ * v.x_ + y_ * v.y_);
	}
	int operator*(const Point& v) const {
		return dot_product(v);
	}
	// for *dot* 
	// http://stackoverflow.com/questions/20861231/designing-a-dot-product-operator-in-c
	// http://codereview.stackexchange.com/questions/23179/named-operators-in-c
	/* TODO: Not working yet
	template<typename LHS>
	struct half_dot {
	LHS lhs;
	half_dot(LHS&& lhs_) :lhs(std::forward<LHS>(lhs_)) {}

	template<typename RHS>
	decltype( dot_product ( std::declval<LHS>(), std::declval<RHS>() ) )
	operator*(RHS&& rhs) const {
	return dot_product(std::forward<LHS>(lhs), std::forward<RHS>(rhs));
	}
	};
	struct dot_t {};
	template<typename LHS>
	half_dot<LHS> operator*(LHS&& lhs, dot_t) {
	return{ std::forward<LHS>(lhs) };
	}
	static dot_t dot;
	*/
	// cos(90) = 0 -> perpendicular lines = dot product 0
	// cos(0) = 1 -> parallel lines = maximum dot product


	int theta(const Point& v) {
		// A dot B = |A||B|Cos(theta)
		return acos(dot_product(v) / (magnitude() + v.magnitude()));
	}

	// cross product of 2 vectors "normally" a vector in +Z direction
	// but for 2D geometry only, considered just a scalar
	// A X B = |A||B|Sin(theta)
	// |theta| = angle between two vectors
	// negative or positive based on "right hand rule"
	// in 2D geometry, if A < 180 degrees clockwise from B, value is positive
	// 1/2 abs(A X B) = area of triangle formed by closing vector A and B
	int cross_product(const Point& v) {
		return (x_ * v.y_ - y_ * v.x_);
	}
	int operator^(const Point& v) {
		return cross_product(v);
	}
	bool operator<(const Point& p) {
		return (x_ < p.x_ && y_ < p.y_);
	}


	// linepoint distance (find shortest distance from point C perpendicular to line formed by AB)
	//     B *
	//      /    
	//   A *--* C
	//  == abs((AB X AC)/|AB|)
	// because area of triangle = base * height / 2
	// area(ABC) =  (AB X AC)/2, base = |AB|
	// height = line distance
	//        = 2 * (AB X AC) / 2 * |AB|
	int linePoint(int aX, int aY, int bX, int bY, int cX, int cY) {
		Point AB(bX - aX, bY - aY);
		Point AC(cX - aX, cY - aY);
		return abs(AB.cross_product(AC) / AB.magnitude());
	}
	// for distance from segment to a point, the shortest path is going to be directly from B to C
	// AB dot BC > 0 implies angle between AB and BC is > -90 < 90, so nearest point is B
	// if BA dot AC > 0, nearest point is A
	// if both dot products are negative, nearest point is somewhere along the segment
	int linePointSegment(int aX, int aY, int bX, int bY, int cX, int cY) {
		Point AB(bX - aX, bY - aY);
		Point BC(cX - bX, cY - bY);
		if (AB.dot_product(BC) > 0) return BC.magnitude();
		Point BA(aX - bX, aY - bY);
		Point AC(cX - aX, cY - aY);
		if (BA.dot_product(AC) > 0) return AC.magnitude();
		return linePoint(aX, aY, bX, bY, cX, cY);
	}

	// point class with operator overloading on points and operator* for dot and ^ for cross 
	double linePointDist(Point A, Point B, Point C, bool isSegment) {
		if (isSegment) {
			int dot1 = (C-B)*(B-A);
			if (dot1 > 0) return (B-C).magnitude();
			int dot2 = (C-A)*(A-B);
			if (dot2 > 0) return (A-C).magnitude();
		}
		return abs((B-A)^(C-A)) / (B-A).magnitude();
	}

	// polygon area (multiple vector points)
	//                *B  
	//              /  \
	//             /    \
	//            /      \
	//           /        \
	//          /          \
	//         /            \ 
	//        /              \ 
	//       *-----------*    \
	//       A        E  |     *C
	//                   |    / 
	//                   |   / 
	//                   |  /
	//                    *
	//                    D
	// Triangulate:
	// ABC + ACD - ADE
	// clockwise + clockwise + counterclockwise
	// AB x AC + AC x AD + AD x AE
	int area(vector<Point> v) {
		int a = 0;
		int N = v.size();
		Point v0 = v[0];
		for (int i = 1; i + 1 < N; i++) {
			Point v1 = v[i] - v0;
			Point v2 = v[i + 1] - v0;
			a += v1 * v2;
		}
		return abs(a / 2);
	}
	// rotate this point around the reference point
	Point rotate(const Point& p, double theta) const {
		// shift our point so p = origin
		Point shifted(x_ - p.x_, y_ - p.y_);
		Point rotated(shifted.x_ * cos(theta) - shifted.y_ * sin(theta), shifted.x_ * sin(theta) + shifted.y_ * cos(theta));
		// shift back
		return Point(rotated.x_ + p.x_, rotated.y_ + p.y_);
	}
	vector<Point> convexHull(vector<Point> &points, bool includeEdge) {
		vector<Point> hull;
		// find left/topmost point index: P
		int len = points.size();
		int P = 0;
		vector<bool> used(len, true);
		for (int i = 1; i < len; i++) {
			if (points[i] < points[P])
				P = i;
		}
		// move 'clockwise', using cross products:
		//  - chose some unused point index: N
		//  - compare to all other point indices: X
		//     - if (X-P) x (N-P) < 0, N = X
		//  - at end of this, now P = N
		int start = P;
		do {
			int N = -1;
			int dist = includeEdge ? INT_MAX : 0;
			for (int i = 0; i < len; i++) {
				// skip previous point
				if (i == P) continue;
				// don't revisit a used point
				if (used[i]) continue;
				// no N yet, set it to i
				if (N == -1) N = i;
				int cross = (points[i] - points[P]) ^ (points[N] - points[P]);
				// d is distance from P to current point
				int d = (points[i] - points[P]) * (points[i] - points[P]);
				if (cross < 0) {
					N = i;
					dist = d;
				}
				else if (cross == 0) {
					// both N and current i are in same direction
					// if includeEdge, pick closest one, otherwise pick furthest
					if (includeEdge && d < dist) {
						dist = d;
						N = i;
					}
					else if (!includeEdge && d > dist) {
						dist = d;
						N = i;
					}
				}
			}
			hull.push_back(points[P]);
		    P = N;
			used[P] = true;
		} while (start != P);
		return hull;
	}
	string pointInPolygon(vector<Point> verts, Point p) {
		int len = verts.size();
		int cnt = 0;
		double x2 = rand() * 1000 + 1000;
		double y2 = rand() * 1000 + 1000;
		for (int i = 0; i < len; i++) {

		}
	}
};
class Line {
public:
	// line: Ax+By=C
	int a_;
	int b_;
	int c_;
	Line(int A, int B, int C) : a_(A), b_(B), c_(C) {}
	Line(Point& p1, Point& p2) {
		a_ = p2.y_ - p1.y_;
		b_ = p1.x_ - p2.x_;
		c_ = a_ * p1.x_ + b_ * p1.y_;
	}
	// 
	Point intersects(const Line& l2) const {
		double det = a_ * l2.b_ - l2.a_ * b_;
		if (det == 0) {
			throw "Parallel";
		}
		Point p(( l2.b_ * c_ - b_ * l2.c_) / det,
			         ( a_ * l2.c_ - l2.a_ * c_) / det);
		return p;
	}
	Point reflection(const Line& refLine, const Point& p) {
		Line perpLine(-refLine.b_, refLine.a_, (-refLine.b_*p.x_ + refLine.a_*p.y_));
		Point inter = refLine.intersects(perpLine);
		return (inter - (p - inter));
	}
};
class LineSegment {
	Point p1_;
	Point p2_;
public:
	LineSegment(Point& p1, Point& p2): p1_(p1), p2_(p2) {}
	// for line segments, (x1, x1) to (x2, y2), get point and validate min (x1, x2) <= x <= max(x1, x2) and same for y
	// due to precision issues with double, need to compare with some tolerance
	Point intersects(const LineSegment& l2) {
		Line line1(p1_, p2_);
		Line line2(l2.p1_, l2.p2_);
		Point p = line1.intersects(line2);
		if (fmin(p1_.x_,p2_.x_) <= p.x_ && 
			p.x_ <= fmax(p1_.x_,p2_.x_)  &&
			fmin(p1_.y_, p2_.y_) <= p.y_ &&
			p.y_ <= fmax(p1_.y_, p2_.y_)  &&
			fmin(l2.p1_.x_, l2.p2_.x_) <= p.x_ &&
			p.x_ <= fmax(l2.p1_.x_, l2.p2_.x_)  &&
			fmin(l2.p1_.y_, l2.p2_.y_) <= p.y_ &&
			p.y_ <= fmax(l2.p1_.y_, l2.p2_.y_)) {
			return p;
		}
		else throw "Don't intersect within segment";
	}
};

class Circle {
	double rad_;
	Point center_;
public:
	// find perpedincular bisectors of XY and YZ
	// then find the point where they intersect, = the center
	Circle(Point& X, Point& Y, Point& Z) {
		Line XY = Line(X, Y);
		Line YZ = Line(Y, Z);
		Point midXY = Point(Y.x_ - ((Y.x_ - X.x_) / 2), Y.y_ - ((Y.y_ - X.y_) / 2));
		Point midYZ = Point(Z.x_ - ((Z.x_ - Y.x_) / 2), Z.y_ - ((Z.y_ - Y.y_) / 2));
		// Ax + By = C
		// perpendicular:
		// -Bx + Ay = D
		Line biSectXY = Line(-XY.b_, XY.a_, (-XY.b_*midXY.x_ + XY.a_*midXY.y_));
		Line biSectYZ = Line(-YZ.b_, YZ.a_, (-YZ.b_*midYZ.x_ + YZ.a_*midYZ.y_));
		center_ = biSectXY.intersects(biSectYZ);
		rad_ = Point(center_.x_ - X.x_, center_.y_ - X.y_).magnitude();
	}
};


int _tmain(int argc, _TCHAR* argv[])
{
	return 0;
}

