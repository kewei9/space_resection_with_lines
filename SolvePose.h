
#ifndef _SOLVE_POSE_FILE_HEADER_
#define _SOLVE_POSE_FILE_HEADER_

#include "powell.h"
#include "SolveDLT.h"

using namespace std;

struct TKBLine
{
	double k;
	double b;
};

void SolvePoseWithLineFeature(TLine2D* linesImg, TPattern pattern, TCamera cam, TPose& pose);


Point2d ProjPerspective(TCamera cam, Point3d pt);

class enginePose : public engine
{
public:
	enginePose();
	~enginePose();

	void SetParameters(TPattern pattern, TCamera cam, TLine2D* imgLines, Mat img);

	TPose GetPoseResult(double* cc);

	virtual double compute_energy(double* current_value);

	void TestPoseValue(TPose pose);

	double GetImageValues(std::vector<Point2d>& points);

private:
	TPattern m_pattern;
	TCamera m_camIn;

	TKBLine m_imgKBLine[6];
	Point3d m_lpt2[6], m_lpt1[6];

	Mat m_distImg;
};


#endif