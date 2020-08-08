
#include "stdafx.h"
#include "SolvePose.h"
#include "gslExMatrix.h"

void GetProjectLinePoints(TCamera cam, TPattern ptn, std::vector<Point2d>& points);

enginePose::enginePose()
{

}

enginePose::~enginePose()
{

}

void enginePose::SetParameters(TPattern pattern, TCamera cam, TLine2D* imgLines, Mat img)
{
	m_pattern = pattern;

	m_lpt2[0] = pattern.ptA;
	m_lpt1[0] = pattern.ptF;

	m_lpt2[1] = pattern.ptC;
	m_lpt1[1] = pattern.ptF;

	m_lpt2[2] = pattern.ptE;
	m_lpt1[2] = pattern.ptF;

	m_lpt2[3] = pattern.ptA;
	m_lpt1[3] = pattern.ptG;

	m_lpt2[4] = pattern.ptB;
	m_lpt1[4] = pattern.ptG;

	m_lpt2[5] = pattern.ptD;
	m_lpt1[5] = pattern.ptG;

	for (int i = 0; i < 6; ++ i)
	{
		double x1 = imgLines[i].pt1.x;
		double y1 = imgLines[i].pt1.y;

		double x2 = imgLines[i].pt2.x;
		double y2 = imgLines[i].pt2.y;

		m_imgKBLine[i].k = (x2 - x1) / (y2 - y1);
		m_imgKBLine[i].b = x1 - m_imgKBLine[i].k * y1;
	}

	convertScaleAbs(img, img, -1.0, 255);

	Mat distImg;
	distanceTransform(img, distImg, CV_DIST_L2, 3);

	double minVal, maxVal;
	minMaxLoc(distImg, &minVal, &maxVal);

	double ra = 255 / maxVal;
	convertScaleAbs(distImg, m_distImg, ra);

	imwrite("D:\\dist_my.bmp", m_distImg);

	m_camIn = cam;
}

TPose enginePose::GetPoseResult(double* cc)
{
	TPose ps = {0};

	double alfa = cc[1];
	double beta = cc[2];
	double gamma= cc[3];

	double mR[9];
	GetCameraRMat(alfa, beta, gamma, mR);

	double matA[36];
	double vecL[12];

	double* mA = matA;
	double k, b, X, Y, Z;

	double x0 = m_camIn.cx;
	double y0 = m_camIn.cy;
	double fx = m_camIn.fx;
	double fy = m_camIn.fy;

	for (int i = 0; i < 6; ++ i)
	{
		k = m_imgKBLine[i].k;
		b = m_imgKBLine[i].b;

		X = m_lpt2[i].x;
		Y = m_lpt2[i].y;
		Z = m_lpt2[i].z;

		double A = mR[0] * X + mR[1] * Y + mR[2] * Z;
		double B = mR[3] * X + mR[4] * Y + mR[5] * Z;
		double C = mR[6] * X + mR[7] * Y + mR[8] * Z;

		mA[0] = fx;
		mA[1] = -k * fy;
		mA[2] = x0 - k * y0 - b;

		vecL[i * 2] = k * y0 * C + k * fy * B + b * C - fx * A - x0 * C; 

		mA += 3;

		X = m_lpt1[i].x;
		Y = m_lpt1[i].y;
		Z = m_lpt1[i].z;

		A = mR[0] * X + mR[1] * Y + mR[2] * Z;
		B = mR[3] * X + mR[4] * Y + mR[5] * Z;
		C = mR[6] * X + mR[7] * Y + mR[8] * Z;

		mA[0] = fx;
		mA[1] = -k * fy;
		mA[2] = x0 - k * y0 - b;

		vecL[i * 2 + 1] =  k * y0 * C + k * fy * B + b * C - fx * A - x0 * C; 

		mA += 3;
	}

	double T[3];
	gslExLeastSquareFit_Linear(matA, vecL, T, 12, 3);

	double t1 = -mR[0] * T[0] - mR[3] * T[1] - mR[6] * T[2];
	double t2 = -mR[1] * T[0] - mR[4] * T[1] - mR[7] * T[2];
	double t3 = -mR[2] * T[0] - mR[5] * T[1] - mR[8] * T[2];

	ps.X0 = t1;
	ps.Y0 = t2;
	ps.Z0 = t3;
	memcpy(ps.matR, mR, sizeof(double) * 9);

	return ps;
}

double enginePose::compute_energy(double* current_value)
{
	double alfa = current_value[1];
	double beta = current_value[2];
	double gamma= current_value[3];

	if (alfa > 180 || beta > 90 || gamma > 90 || alfa < 0 || beta < -90 || gamma < -90)
		return 100000000000;

	double mR[9];
	GetCameraRMat(alfa, beta, gamma, mR);

	double matA[36];
	double vecL[12];

	double* mA = matA;
	double k, b, X, Y, Z;

	double x0 = m_camIn.cx;
	double y0 = m_camIn.cy;
	double fx = m_camIn.fx;
	double fy = m_camIn.fy;

	for (int i = 0; i < 6; ++ i)
	{
		k = m_imgKBLine[i].k;
		b = m_imgKBLine[i].b;

		X = m_lpt2[i].x;
		Y = m_lpt2[i].y;
		Z = m_lpt2[i].z;

		double A = mR[0] * X + mR[1] * Y + mR[2] * Z;
		double B = mR[3] * X + mR[4] * Y + mR[5] * Z;
		double C = mR[6] * X + mR[7] * Y + mR[8] * Z;

		mA[0] = fx;
		mA[1] = -k * fy;
		mA[2] = x0 - k * y0 - b;

		vecL[i * 2] = k * y0 * C + k * fy * B + b * C - fx * A - x0 * C; 

		mA += 3;

		X = m_lpt1[i].x;
		Y = m_lpt1[i].y;
		Z = m_lpt1[i].z;

		A = mR[0] * X + mR[1] * Y + mR[2] * Z;
		B = mR[3] * X + mR[4] * Y + mR[5] * Z;
		C = mR[6] * X + mR[7] * Y + mR[8] * Z;

		mA[0] = fx;
		mA[1] = -k * fy;
		mA[2] = x0 - k * y0 - b;

		vecL[i * 2 + 1] =  k * y0 * C + k * fy * B + b * C - fx * A - x0 * C; 

		mA += 3;
	}

	double T[3];
	gslExLeastSquareFit_Linear(matA, vecL, T, 12, 3);

	double t1 = -mR[0] * T[0] - mR[3] * T[1] - mR[6] * T[2];
	double t2 = -mR[1] * T[0] - mR[4] * T[1] - mR[7] * T[2];
	double t3 = -mR[2] * T[0] - mR[5] * T[1] - mR[8] * T[2];

	m_camIn.camPose.X0 = t1;
	m_camIn.camPose.Y0 = t2;
	m_camIn.camPose.Z0 = t3;
	memcpy(m_camIn.camPose.matR, mR, sizeof(double) * 9);

	std::vector<Point2d> points;
	GetProjectLinePoints(m_camIn, m_pattern, points);

	if (points.size() < 1200 * 6)
		return 10000000000;

	double energy = GetImageValues(points);

	return energy;
}

void enginePose::TestPoseValue(TPose pose)
{

}

double enginePose::GetImageValues(std::vector<Point2d>& points)
{
	uchar* dat = (uchar*)(m_distImg.data);

	int w = m_distImg.cols;

	double en = 0.0f;
	for (int i = 0; i < points.size(); ++ i)
	{
		en += *(dat + int(points[i].y + 0.5f) * w + int(points[i].x + 0.5f));
	}

// 	Mat imgTemp = m_distImg.clone();
// 	Mat clrTemp;
// 	cvtColor(imgTemp, clrTemp, CV_GRAY2RGB);
// 
// 	for (int i = 0; i < points.size(); ++ i)
// 	{
// 		Vec3b &p = clrTemp.at<Vec3b>(int(points[i].y + 0.5f), int(points[i].x + 0.5f));
// 		p[0] = 0;
// 		p[1] = 0;
// 		p[2] = 255;
// 	}
// 
// 	imwrite("D:\\init_val_img.bmp", clrTemp);

	return en;
}

Point2d ProjPerspective(TCamera cam, Point3d pt)
{
	Point3d p1 = pt;

	double fx = cam.fx;
	double fy = cam.fy;

	double cx = cam.cx;
	double cy = cam.cy;

	double Xs = cam.camPose.X0;
	double Ys = cam.camPose.Y0;
	double Zs = cam.camPose.Z0;

	double* mR = cam.camPose.matR;

	double X = mR[0] * (p1.x - Xs) + mR[1] * (p1.y - Ys) + mR[2] * (p1.z - Zs);
	double Y = mR[3] * (p1.x - Xs) + mR[4] * (p1.y - Ys) + mR[5] * (p1.z - Zs);
	double Z = mR[6] * (p1.x - Xs) + mR[7] * (p1.y - Ys) + mR[8] * (p1.z - Zs);

	double x1 = cx + fx * X / Z;
	double y1 = cy + fy * Y / Z;

	Point2d res;
	res.x = x1;
	res.y = y1;

	return res;
}

void GetLinePoints(TCamera cam, TLine3D lineSpace, std::vector<Point2d>& points)
{
	Point2d pt1 = ProjPerspective(cam, lineSpace.pt1);
	Point2d pt2 = ProjPerspective(cam, lineSpace.pt2);

	double x1 = pt1.x;
	double y1 = pt1.y;
	double x2 = pt2.x;
	double y2 = pt2.y;

	int w = 1280;
	int h = 1280;

	int nk = abs(x2 - x1);
	if (abs(y2 - y1) > nk)
		nk = abs(y2 - y1);

	float x_incre = (x2 - x1) * 1.0f / nk;
	float y_incre = (y2 - y1) * 1.0f / nk;

	float x = (float)(x1);
	float y = (float)(y1);

	for (int i = 0; i < nk; ++ i)
	{	
		if (x > w || y > h || x < 0 || y < 0)
		{
			x += x_incre;
			y += y_incre;
			continue;
		}

		x += x_incre;
		y += y_incre;

		Point2d pt(x, y);
		points.push_back(pt);
	}
}

void GetProjectLinePoints(TCamera cam, TPattern ptn, std::vector<Point2d>& points)
{
	points.clear();

	TLine3D spcLineTmp;

	spcLineTmp.pt1 = ptn.ptA;
	spcLineTmp.pt2 = ptn.ptF;
	GetLinePoints(cam, spcLineTmp, points);

	spcLineTmp.pt1 = ptn.ptC;
	spcLineTmp.pt2 = ptn.ptF;
	GetLinePoints(cam, spcLineTmp, points);

	spcLineTmp.pt1 = ptn.ptE;
	spcLineTmp.pt2 = ptn.ptF;
	GetLinePoints(cam, spcLineTmp, points);

	spcLineTmp.pt1 = ptn.ptA;
	spcLineTmp.pt2 = ptn.ptG;
	GetLinePoints(cam, spcLineTmp, points);

	spcLineTmp.pt1 = ptn.ptB;
	spcLineTmp.pt2 = ptn.ptG;
	GetLinePoints(cam, spcLineTmp, points);

	spcLineTmp.pt1 = ptn.ptD;
	spcLineTmp.pt2 = ptn.ptG;
	GetLinePoints(cam, spcLineTmp, points);
}


void multiplyRMat(double* R1, double* R2, double* R)
{
	R[0] = R1[0] * R2[0] + R1[1] * R2[3] + R1[2] * R2[6];
	R[1] = R1[0] * R2[1] + R1[1] * R2[4] + R1[2] * R2[7];
	R[2] = R1[0] * R2[2] + R1[1] * R2[5] + R1[2] * R2[8];

	R[3] = R1[3] * R2[0] + R1[4] * R2[3] + R1[5] * R2[6];
	R[4] = R1[3] * R2[1] + R1[4] * R2[4] + R1[5] * R2[7];
	R[5] = R1[3] * R2[2] + R1[4] * R2[5] + R1[5] * R2[8];

	R[6] = R1[6] * R2[0] + R1[7] * R2[3] + R1[8] * R2[6];
	R[7] = R1[6] * R2[1] + R1[7] * R2[4] + R1[8] * R2[7];
	R[8] = R1[6] * R2[2] + R1[7] * R2[5] + R1[8] * R2[8];
}


void GetCameraRMat(double alfa, double beta, double sita, double* matR)
{
	double sin_alfa = sin(alfa / 180 * CV_PI);
	double cos_alfa = cos(alfa / 180 * CV_PI);

	double sin_beta = sin(beta / 180 * CV_PI);
	double cos_beta = cos(beta / 180 * CV_PI);

	double sin_sita = sin(sita / 180 * CV_PI);
	double cos_sita = cos(sita / 180 * CV_PI);

	double m1[] = {1.0, 0.0, 0.0,  0.0, cos_alfa, sin_alfa, 0.0, -sin_alfa, cos_alfa};
	double m2[] = {cos_beta, 0.0, sin_beta, 0.0, 1.0, 0.0, -sin_beta, 0.0, cos_beta};
	double m3[] = {cos_sita, sin_sita, 0.0, -sin_sita, cos_sita, 0.0, 0.0, 0.0, 1.0};

	double matT[9];
	multiplyRMat(m1, m2, matT);
	multiplyRMat(matT, m3, matR);
}

void SolvePoseWithLineFeature(TLine2D* linesImg, TPattern pattern, TCamera cam, TPose& pose)
{
	double c[4] = {0, 180.0f, 0.0f, 8.0f};

	Mat image = imread("D:\\org_image.bmp", CV_LOAD_IMAGE_GRAYSCALE);

	enginePose myEngine;
	myEngine.SetParameters(pattern, cam, linesImg, image);

	myEngine.register_engine();

	optimization_powell powell;
	double res = powell.optimize(c, 3, 0.00000000000001);

	pose = myEngine.GetPoseResult(c);

}