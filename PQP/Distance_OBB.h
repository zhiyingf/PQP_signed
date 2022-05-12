#pragma once
#include "..\Model3D\BaseModel.h"
#include "..\Model3D\Polygon2Model3D.hpp"

#include <iostream>
#include <Eigen/dense>
#include <algorithm>
#include <CGAL/Point_2.h>
#include <CGAL/Point_3.h>
#include <CGAL/Polygon_2.h>
#include <CGAL\Exact_predicates_inexact_constructions_kernel.h>

using namespace std;

#include "PQP.h"

namespace PQP
{
	typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
	typedef CGAL::Polygon_2<K> Polygon_2;
	typedef CGAL::Point_2<K> Point_2;
	typedef CGAL::Point_3<K> Point_3;
	typedef CGAL::Vector_3<K> Vector_3;

	using namespace Model3D;
	class Distance_OBB
	{
		CBaseModel m_model;
		PQP_Model m_pqp_model;
	public:
		struct QueryResult
		{
			int pos_flag;
			// pos_flag == 0:  inside
			// pos_flag == 1: edge
			// pos_flag == 2: edge
			// pos_flag == 3: edge
			// pos_flag == 4: vertex
			// pos_flag == 5: vertex
			// pos_flag == 6: vertex		

			int triangleID;
			double distance;
			CPoint3D closestPnt;
			QueryResult() {}
			QueryResult(int triangleID, double distance, const CPoint3D &closestPnt, int pos_flag)
				: triangleID(triangleID), distance(distance), closestPnt(closestPnt), pos_flag(pos_flag)
			{
			}
			bool operator<(const QueryResult& other) const
			{
				if (distance < other.distance)
					return true;
				if (distance > other.distance)
					return false;
				if (triangleID < other.triangleID)
					return true;
				if (triangleID > other.triangleID)
					return false;
				return false;
			}
		};
	public:
		Distance_OBB();
		Distance_OBB(CBaseModel& model);
		QueryResult Query(const CPoint3D& pt);
		vector<QueryResult> Query(const CPoint3D& pt, const CPoint3D& dir);
		vector<QueryResult> SegQuery(const CPoint3D& source, const CPoint3D& target);
		bool LyingInside(const CPoint3D& pt);
		double SignedDistance(const CPoint3D& pt);
		double SignedDistance(const CPoint3D& pt, CPoint3D& gradient);

		template<typename K>
		double SignedDistance(const CGAL::Point_2<K>& pt)
		{
			return SignedDistance(CPoint3D(pt.x(), pt.y(), 0));
		}
		template<typename K>
		double SignedDistance(const CGAL::Point_2<K>& pt, CGAL::Point_2<K>& gradient)
		{
			CPoint3D gradient3d;
			auto distance = SignedDistance(CPoint3D(pt.x(), pt.y(), 0), gradient3d);
			gradient = Point_2(gradient3d.x, gradient3d.y);
			return distance;
		}
		template<typename K>
		double SignedDistance(const CGAL::Point_3<K>& pt)
		{
			return SignedDistance(CPoint3D(pt.x(), pt.y(), pt.z()));
		}
		template<typename K>
		double SignedDistance(const CGAL::Point_3<K>& pt, CGAL::Point_3<K>& gradient)
		{
			CPoint3D gradient3d;
			auto distance = SignedDistance(CPoint3D(pt.x(), pt.y(), pt.z()), gradient3d);
			gradient = Point_3(gradient3d.x, gradient3d.y, gradient3d.z);
			return distance;
		}
	};
}