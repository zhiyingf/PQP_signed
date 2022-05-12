#include "Distance_OBB.h"

namespace PQP
{
	/*typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
	typedef CGAL::Polygon_2<K> Polygon_2;
	typedef CGAL::Point_2<K> Point_2;
	typedef CGAL::Point_3<K> Point_3;
	typedef CGAL::Vector_3<K> Vector_3;

	using namespace Model3D;*/

	Distance_OBB::Distance_OBB()
	{	}

	Distance_OBB::Distance_OBB(CBaseModel& model)
		: m_model(model)
	{
		m_pqp_model.BeginModel();
		PQP_REAL p1[3], p2[3], p3[3];
		for (int i = 0; i < m_model.GetNumOfFaces(); ++i)
		{
			int vid = m_model.Face(i)[0];
			p1[0] = (PQP_REAL)(m_model.Vert(vid).x);
			p1[1] = (PQP_REAL)(m_model.Vert(vid).y);
			p1[2] = (PQP_REAL)(m_model.Vert(vid).z);
			vid = m_model.Face(i)[1];
			p2[0] = (PQP_REAL)(m_model.Vert(vid).x);
			p2[1] = (PQP_REAL)(m_model.Vert(vid).y);
			p2[2] = (PQP_REAL)(m_model.Vert(vid).z);
			vid = m_model.Face(i)[2];
			p3[0] = (PQP_REAL)(m_model.Vert(vid).x);
			p3[1] = (PQP_REAL)(m_model.Vert(vid).y);
			p3[2] = (PQP_REAL)(m_model.Vert(vid).z);
			m_pqp_model.AddTri(p1, p2, p3, i);
		}
		m_pqp_model.EndModel();
	}

	bool Distance_OBB::LyingInside(const CPoint3D& pt)
	{
		auto queryRes = Query(pt);
		double dis1 = (queryRes.closestPnt - m_model.Vert(m_model.Face(queryRes.triangleID)[0])).Len();
		double dis2 = (queryRes.closestPnt - m_model.Vert(m_model.Face(queryRes.triangleID)[1])).Len();
		double dis3 = (queryRes.closestPnt - m_model.Vert(m_model.Face(queryRes.triangleID)[2])).Len();
		CPoint3D averageNormal(0, 0, 0);
		averageNormal = averageNormal + m_model.Normal(m_model.Face(queryRes.triangleID)[0]) * (dis2 + dis3);
		averageNormal = averageNormal + m_model.Normal(m_model.Face(queryRes.triangleID)[1]) * (dis1 + dis3);
		averageNormal = averageNormal + m_model.Normal(m_model.Face(queryRes.triangleID)[2]) * (dis1 + dis2);
		averageNormal.Normalize();
		auto delta = pt - queryRes.closestPnt;
		delta.Normalize();
		if ((delta ^ averageNormal) >= 0)
			return false;
		return true;
	}

	Distance_OBB::QueryResult Distance_OBB::Query(const CPoint3D& pt)
	{
		PQP_DistanceResult dres;
		dres.last_tri = m_pqp_model.last_tri;
		PQP_REAL p[3];
		p[0] = pt.x;	p[1] = pt.y;	p[2] = pt.z;
		PQP_Distance(&dres, &m_pqp_model, p, 0.0, 0.0);
		Distance_OBB::QueryResult qresult;
		qresult.closestPnt.x = dres.p1[0];
		qresult.closestPnt.y = dres.p1[1];
		qresult.closestPnt.z = dres.p1[2];
		qresult.triangleID = dres.last_tri->id;
		qresult.distance = dres.distance;
		qresult.pos_flag = dres.pos_flag;

		//modify to signed distance
		/*double dis1 = (qresult.closestPnt - m_model.Vert(m_model.Face(qresult.triangleID)[0])).Len();
		double dis2 = (qresult.closestPnt - m_model.Vert(m_model.Face(qresult.triangleID)[1])).Len();
		double dis3 = (qresult.closestPnt - m_model.Vert(m_model.Face(qresult.triangleID)[2])).Len();
		CPoint3D averageNormal(0, 0, 0);
		averageNormal = averageNormal + m_model.Normal(m_model.Face(qresult.triangleID)[0]) * (dis2 + dis3);
		averageNormal = averageNormal + m_model.Normal(m_model.Face(qresult.triangleID)[1]) * (dis1 + dis3);
		averageNormal = averageNormal + m_model.Normal(m_model.Face(qresult.triangleID)[2]) * (dis1 + dis2);
		averageNormal.Normalize();
		auto delta = pt - qresult.closestPnt;
		delta.Normalize();
		if ((delta ^ averageNormal) >= 0)
			qresult.distance = qresult.distance;
		else
			qresult.distance = -qresult.distance;*/
		//end

		return qresult;
	}

	pair<CPoint3D, double> GetIntersectionPoint(const CPoint3D& pt1, const CPoint3D& pt2, const CPoint3D& pt3,
		const CPoint3D& pt, const CPoint3D& dir)
	{
		Eigen::MatrixXd M(4, 4);
		M << pt1.x, pt2.x, pt3.x, -dir.x,
			pt1.y, pt2.y, pt3.y, -dir.y,
			pt1.z, pt2.z, pt3.z, -dir.z,
			1, 1, 1, 0;
		Eigen::VectorXd right(4);
		right << pt.x, pt.y, pt.z, 1;
		Eigen::VectorXd res = M.inverse() * right;
		double t = res(3);
		return make_pair(pt + t * dir, t);
	}
	vector<Distance_OBB::QueryResult> Distance_OBB::Query(const CPoint3D& pt, const CPoint3D& dir)
	{
		double M = 1000;
		CPoint3D start = pt + 1e-4 * dir;
		CPoint3D end = pt + M * dir;
		PQP_Model singleTriangle;
		singleTriangle.BeginModel();
		PQP_REAL p1[3], p2[3], p3[3];
		p1[0] = (PQP_REAL)(start.x);
		p1[1] = (PQP_REAL)(start.y);
		p1[2] = (PQP_REAL)(start.z);
		p2[0] = (PQP_REAL)(end.x);
		p2[1] = (PQP_REAL)(end.y);
		p2[2] = (PQP_REAL)(end.z);
		p3[0] = (PQP_REAL)(end.x);
		p3[1] = (PQP_REAL)(end.y);
		p3[2] = (PQP_REAL)(end.z);
		singleTriangle.AddTri(p1, p2, p3, 0);
		singleTriangle.EndModel();
		PQP_CollideResult result;
		PQP_REAL R1[3][3] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
		PQP_REAL T1[3] = { 0, 0, 0 };
		PQP_REAL R2[3][3] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
		PQP_REAL T2[3] = { 0, 0, 0 };
		int ret = PQP_Collide(&result, R1, T1, &m_pqp_model, R2, T2, &singleTriangle);
		vector<Distance_OBB::QueryResult> intersections;
		for (int i = 0; i < result.num_pairs; ++i)
		{
			pair<CPoint3D, double> intersection = GetIntersectionPoint(m_model.Vert(m_model.Face(result.Id1(i))[0]),
				m_model.Vert(m_model.Face(result.Id1(i))[1]),
				m_model.Vert(m_model.Face(result.Id1(i))[2]),
				start, dir);
			intersections.push_back(Distance_OBB::QueryResult(result.Id1(i),
				intersection.second, intersection.first, -1));
		}
		std::sort(intersections.begin(), intersections.end());
		return intersections;
	}
	pair<CPoint3D, double> GetIntersectionPoint2(const CPoint3D& pt1, const CPoint3D& pt2, const CPoint3D& pt3,
		const CPoint3D& pt, const CPoint3D& ptt)
	{
		CPoint3D dir = ptt - pt;
		dir.Normalize();

		Eigen::MatrixXd M(4, 4);
		M << pt1.x, pt2.x, pt3.x, -dir.x,
			pt1.y, pt2.y, pt3.y, -dir.y,
			pt1.z, pt2.z, pt3.z, -dir.z,
			1, 1, 1, 0;
		Eigen::VectorXd right(4);
		right << pt.x, pt.y, pt.z, 1;
		Eigen::VectorXd res = M.inverse() * right;
		double t = res(3);
		return make_pair(pt + t * dir, t);
	}
	vector<Distance_OBB::QueryResult> Distance_OBB::SegQuery(const CPoint3D & source, const CPoint3D & target)
	{
		CPoint3D start = source;
		CPoint3D end = target;
		PQP_Model singleTriangle;
		singleTriangle.BeginModel();
		PQP_REAL p1[3], p2[3], p3[3];
		p1[0] = (PQP_REAL)(start.x);
		p1[1] = (PQP_REAL)(start.y);
		p1[2] = (PQP_REAL)(start.z);
		p2[0] = (PQP_REAL)(end.x);
		p2[1] = (PQP_REAL)(end.y);
		p2[2] = (PQP_REAL)(end.z);
		p3[0] = (PQP_REAL)(end.x);
		p3[1] = (PQP_REAL)(end.y);
		p3[2] = (PQP_REAL)(end.z);
		singleTriangle.AddTri(p1, p2, p3, 0);
		singleTriangle.EndModel();
		PQP_CollideResult result;
		PQP_REAL R1[3][3] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
		PQP_REAL T1[3] = { 0, 0, 0 };
		PQP_REAL R2[3][3] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
		PQP_REAL T2[3] = { 0, 0, 0 };
		int ret = PQP_Collide(&result, R1, T1, &m_pqp_model, R2, T2, &singleTriangle);
		vector<Distance_OBB::QueryResult> intersections;
		for (int i = 0; i < result.num_pairs; ++i)
		{
			pair<CPoint3D, double> intersection = GetIntersectionPoint2(m_model.Vert(m_model.Face(result.Id1(i))[0]),
				m_model.Vert(m_model.Face(result.Id1(i))[1]),
				m_model.Vert(m_model.Face(result.Id1(i))[2]),
				start, end);
			intersections.push_back(Distance_OBB::QueryResult(result.Id1(i),
				intersection.second, intersection.first, -1));
		}
		std::sort(intersections.begin(), intersections.end());
		return intersections;
	}

	double Distance_OBB::SignedDistance(const CPoint3D& pt)
	{
		//inside <= 0; outside >= 0
		auto queryRes = Query(pt);
		double dis1 = (queryRes.closestPnt - m_model.Vert(m_model.Face(queryRes.triangleID)[0])).Len();
		double dis2 = (queryRes.closestPnt - m_model.Vert(m_model.Face(queryRes.triangleID)[1])).Len();
		double dis3 = (queryRes.closestPnt - m_model.Vert(m_model.Face(queryRes.triangleID)[2])).Len();
		CPoint3D averageNormal(0, 0, 0);
		averageNormal = averageNormal + m_model.Normal(m_model.Face(queryRes.triangleID)[0]) * (dis2 + dis3);
		averageNormal = averageNormal + m_model.Normal(m_model.Face(queryRes.triangleID)[1]) * (dis1 + dis3);
		averageNormal = averageNormal + m_model.Normal(m_model.Face(queryRes.triangleID)[2]) * (dis1 + dis2);
		averageNormal.Normalize();
		auto delta = pt - queryRes.closestPnt;
		delta.Normalize();
		if ((delta ^ averageNormal) >= 0)
			return queryRes.distance;
		return -queryRes.distance;
	}

	double Distance_OBB::SignedDistance(const CPoint3D& pt, CPoint3D& gradient)
	{
		//inside <= 0; outside >= 0
		auto queryRes = Query(pt);
		double dis1 = (queryRes.closestPnt - m_model.Vert(m_model.Face(queryRes.triangleID)[0])).Len();
		double dis2 = (queryRes.closestPnt - m_model.Vert(m_model.Face(queryRes.triangleID)[1])).Len();
		double dis3 = (queryRes.closestPnt - m_model.Vert(m_model.Face(queryRes.triangleID)[2])).Len();
		CPoint3D averageNormal(0, 0, 0);
		averageNormal = averageNormal + m_model.Normal(m_model.Face(queryRes.triangleID)[0]) * (dis2 + dis3);
		averageNormal = averageNormal + m_model.Normal(m_model.Face(queryRes.triangleID)[1]) * (dis1 + dis3);
		averageNormal = averageNormal + m_model.Normal(m_model.Face(queryRes.triangleID)[2]) * (dis1 + dis2);
		averageNormal.Normalize();
		auto delta = pt - queryRes.closestPnt;
		delta.Normalize();
		if ((delta ^ averageNormal) >= 0)//outside
		{
			gradient = pt - queryRes.closestPnt;
			gradient.Normalize();
			return queryRes.distance;
		}
		gradient = queryRes.closestPnt - pt;
		gradient.Normalize();
		return -queryRes.distance;//inside
	}
}