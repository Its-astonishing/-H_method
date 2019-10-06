#include "CH.h"


double FTEST(double x, double v)
{
	return v;
}

double FBASIC(double x, double v, double a, double s)
{
	const double g = 9.8;
	const double pi = 3.14159265;
	return(-0.6 * s * sqrt(2 * g) * (1 / pow(tan(0.5 * a), 2)) * (1 / pi) * (1 / pow(sqrt(v), 3)));
}

double RK4_Test( double Xn, double Vn, double hn)
{
	double V; // Vn+1
	double k1, k2, k3, k4;
    
	k1 = FTEST(Xn, Vn);
	k2 = FTEST(Xn + (hn / 2), Vn + (hn / 2) * k1);
	k3 = FTEST(Xn + (hn / 2), Vn + (hn / 2) * k2);
	k4 = FTEST(Xn + hn, Vn + hn * k3);
	
	V = Vn + (hn / 6) * (k1 + (2 * k2) + (2 * k3) + k4);
	
	return V;
}

double RK4_BASIC(double Xn, double Vn, double hn, double a, double s)
{
	double V; // Vn+1
	double k1, k2, k3, k4;

	k1 = FBASIC(Xn, Vn, a, s);
	k2 = FBASIC(Xn + (hn / 2), Vn + (hn / 2) * k1,a,s);
	k3 = FBASIC(Xn + (hn / 2), Vn + (hn / 2) * k2, a, s);
	k4 = FBASIC(Xn + hn, Vn + hn * k3, a, s);

	V = Vn + (hn / 6) * (k1 + (2 * k2) + (2 * k3) + k4);

	return V;
}

std::vector<Point> StartMethodTest(int n, Point p0)
{
	std::vector<Point> points;
	points.push_back(p0);
	double h = p0.h;
	int i = 1;
	double xtmp;
	Point p;
	p.er=p0.er;
	double c = p0.vn - exp(p0.xn);
	while (i < n + 1)
	{
		p.num = i;
		p.h = h;
		p.xn = points[i - 1].xn + h;
		p.vn = RK4_Test(points[i - 1].xn, points[i - 1].vn, h);
		xtmp = p.xn - h / 2;
		p.v15 = RK4_Test(points[i - 1].xn, points[i - 1].vn, h / 2);
		p.v2 = RK4_Test(xtmp, p.v15, h / 2);
		p.S = (p.vn - p.v2) / 7;
		if (abs(p.S) > points[i - 1].er)
		{
			h = h / 2;
			continue;
		}
		if (abs(p.S) < points[i - 1].er / 32)
		{
			h = 2*h;
			p.u = exp(p.xn)+c;
			p.Gerror = p.u - p.vn;
			points.push_back(p);
			i++;
			continue;
		}
		p.u = exp(p.xn)+c;
		p.Gerror = p.u - p.vn;
		points.push_back(p);
		i++;
	}
	return points;
}

std::vector<Point> StartMethodBasic(int n, Point p0, double a, double s)
{
	std::vector<Point> points;
	points.push_back(p0);
	double h = p0.h;
	int i = 1;
	double xtmp;
	Point p;
	p.er = p0.er;
	while (i < n + 1)
	{
		p.num = i;
		p.h = h;
		p.xn = points[i - 1].xn + h;
		p.vn = RK4_BASIC(points[i - 1].xn, points[i - 1].vn, h,a,s);
		xtmp = p.xn - h / 2;
		p.v15 = RK4_BASIC(points[i - 1].xn, points[i - 1].vn, h / 2,a,s);
		p.v2 = RK4_BASIC(xtmp, p.v15, h / 2,a,s);
		p.S = (p.vn - p.v2) / 7;
		if (abs(p.S) > points[i - 1].er)
		{
			h = h / 2;
			continue;
		}
		if (abs(p.S) < points[i - 1].er / 32)
		{
			h = 2 * h;
			points.push_back(p);
			i++;
			continue;
		}
		points.push_back(p);
		i++;
	}
	return points;
}


