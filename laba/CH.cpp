#include <iostream>
#include <cmath>
#include <vector>
#include <functional>

struct PPaint
{
	int num;
	double u;
	double xn, vn; //
	double Gerror; //глобальная погрешность
	double v2, v15; // v удвоенное и v+1/2
	double S, er; // Контрольное слагаемое и оценка лин погрешности
	double h;
	int sub; 
	PPaint()
	{
		num = 0;
		sub = u = xn = vn = Gerror = v2 = v15 = S = er = 0;
		h = 0.1;
	}
};

struct PPaint1
{
	int num;
	double xn, v1, v2; //
	double v22, v21, v115, v215;
	double S, er;
	double h;
	int sub;
	PPaint1()
	{
		num = 0;
		sub = xn = v1 = v2 = v21 = v22 = v115 = v215 = S=er = 0;
		h = 0.1;
	}
};

struct boundaryProblemPoint
{
	int number = 0;
	double step = 0;
	double x = 0;
	double numeric_sol_value			 = 0;
	double numeric_sol_value_doubled	 = 0;
	double numeric_sol_value_halfstepped = 0;
	double global_error			   = 0; 
	double control_value		   = 0;
	double linear_error_evaluation = 0;
	double a = 0;
	double phi = 0;
	double d = 0;
};


inline double max(double x1, double x2)
{
	if (x1 > x2)
		return x1;
	return x2;
}

inline double FTEST(double x, double v)
{
	return v;
}

inline double FBASIC1(double x, double v)
{
	return (x/(1+x*x))*v*v+v-v*v*v*sin(10*x);
}

inline double FBASIC_2_1(double x, double v1, double v2)
{
	return v2;
}

inline double FBASIC_2_2(double x, double v1, double v2, double a,double b)
{
	return a*v2+b*sin(v1);
}

inline double RK4_1(double x, double v1, double v2, double h)
{
	double V; // Vn+1
	double k1, k2, k3, k4;

	k1 = FBASIC_2_1(x, v1, v2);
	k2 = FBASIC_2_1(x + (h / 2.0), v1 + (h / 2.0) * k1, v2 + (h / 2.0) * k1);
	k3 = FBASIC_2_1(x + (h / 2.0), v1 + (h / 2.0) * k2, v2 + (h / 2.0) * k2);
	k4 = FBASIC_2_1(x + h, v1 + h * k3, v2 + h * k3);

	V = v1 + (h / 6.0) * (k1 + (2.0 * k2) + (2.0 * k3) + k4);

	return V;
}

inline double RK4_2(double x, double v1, double v2, double h, double a, double b)
{
	double V; // Vn+1
	double k1, k2, k3, k4;

	k1 = FBASIC_2_2(x, v1, v2, a, b);
	k2 = FBASIC_2_2(x + (h / 2.0), v1 + (h / 2.0) * k1, v2 + (h / 2.0) * k1, a, b);
	k3 = FBASIC_2_2(x + (h / 2.0), v1 + (h / 2.0) * k2, v2 + (h / 2.0) * k2, a, b);
	k4 = FBASIC_2_2(x + h, v1 + h * k3, v2 + h * k3, a, b);

	V = v2 + (h / 6.0) * (k1 + (2.0 * k2) + (2.0 * k3) + k4);

	return V;
}


inline double RK4_Test(double Xn, double Vn, double hn)
{
	double V; // Vn+1
	double k1, k2, k3, k4;

	k1 = FTEST(Xn, Vn);
	k2 = FTEST(Xn + (hn / 2.0), Vn + (hn / 2.0) * k1);
	k3 = FTEST(Xn + (hn / 2.0), Vn + (hn / 2.0) * k2);
	k4 = FTEST(Xn + hn, Vn + hn * k3);

	V = Vn + (hn / 6.0) * (k1 + (2.0 * k2) + (2.0 * k3) + k4);

	return V;
}

inline double RK4_BASIC1(double Xn, double Vn, double hn)
{
	double V; // Vn+1
	double k1, k2, k3, k4;

	k1 = FBASIC1(Xn, Vn);
	k2 = FBASIC1(Xn + (hn / 2.0), Vn + (hn / 2.0) * k1);
	k3 = FBASIC1(Xn + (hn / 2.0), Vn + (hn / 2.0) * k2);
	k4 = FBASIC1(Xn + hn, Vn + hn * k3);

	V = Vn + (hn / 6.0) * (k1 + (2.0 * k2) + (2.0 * k3) + k4);

	return V;
}

inline std::vector<PPaint> StartMethodTest(int n, PPaint p0,double control, double range, double xmax)
{
	std::vector<PPaint> points;
	p0.u = p0.vn;
	points.push_back(p0);
	double h = p0.h;
	int i = 1;
	double xtmp;
	PPaint p;
	double c = p0.vn/exp(p0.xn);
	while (i < n + 1)
	{
		if (points[i - 1].xn > xmax - range)
			break;
		p.num = i;
		p.h = h;
		p.xn = points[i - 1].xn + h;
		if (p.xn > xmax)
		{
			h = h / 2;
			continue;
		}
		p.vn = RK4_Test(points[i - 1].xn, points[i - 1].vn, h);
		xtmp = p.xn - h / 2;
		p.v15 = RK4_Test(points[i - 1].xn, points[i - 1].vn, h / 2.0);
		p.v2 = RK4_Test(xtmp, p.v15, h / 2.0);
		p.S = (p.vn - p.v2) / 15.0;
		if (abs(p.S) > control)
		{
			p.sub = -2;
			h = h / 2.0;
			continue;
		}
		if (abs(p.S) < control/ 32.0)
		{
			p.sub = 2;
			h = 2.0 * h;
			p.u = c*exp(p.xn);
			p.Gerror = p.u - p.vn;
			points.push_back(p);
			i++;
			continue;
		}
		p.u = c*exp(p.xn);
		p.Gerror = p.u - p.vn;
		points.push_back(p);
		i++;
	}
	return points;
}

inline std::vector<PPaint> Test(int n, PPaint p0, double control, double range, double xmax)
{
	std::vector<PPaint> points;
	p0.u = p0.vn;
	points.push_back(p0);
	double h = p0.h;
	int i = 1;
	PPaint p;
	double c = p0.vn/exp(p0.xn);
	while (i < n + 1)
	{
		if (points[i - 1].xn > xmax - range)
			break;
		p.num = i;
		p.h = h;
		p.xn = points[i - 1].xn + h;
		if (p.xn > xmax)
		{
			h = h / 2;
			continue;
		}
		p.vn = RK4_Test(points[i - 1].xn, points[i - 1].vn, h);
		p.u = c*exp(p.xn);
		p.Gerror = p.u - p.vn;
		points.push_back(p);
		i++;
	}
	return points;
}

inline std::vector<PPaint> Basic(int n, PPaint p0, double control, double range, double xmax)
{
	std::vector<PPaint> points;
	points.push_back(p0);
	double h = p0.h;
	int i = 1;
	PPaint p;
	while (i < n + 1)
	{
		if (points[i - 1].xn > xmax - range)
			break;
		p.num = i;
		p.h = h;
		p.xn = points[i - 1].xn + h;
		if (p.xn > xmax)
		{
			h = h / 2;
			continue;
		}
		p.vn = RK4_BASIC1(points[i - 1].xn, points[i - 1].vn, h);
		points.push_back(p);
		i++;
	}
	return points;
}

template<typename Function>
static inline double calc_integral(double &lower_bound, double &upper_bound, Function function_to_be_integrated)
{
	double result_value = 0;

	return result_value;
}

template<typename Function>
static inline double calc_a_value(double xi0, double xi1, Function k_function)
{
	double a_value = 0;
	return return a_value;
}

template<typename Function>
static inline double calc_phi_value(double xi0, double xi1, Function f_function)
{
	double phi_value = 0;
	return return phi_value;
}

template<typename Function>
static inline double calc_d_value(double xi0, double xi1, Function q_function)
{
	double d_value = 0;
	return return d_value;
}

inline std::vector<boundaryProblemPoint> StartMethodBasic1(int n, boundaryProblemPoint p0, double control, double range, double xmax, double gap_point)
{
	auto k1 = [](double &x_value) -> double
	{
		return sqrt(x_value);
	};
	auto k2 = [](double &x_value) -> double
	{
		return x_value + 1;
	};
	auto f1 = [](double &x_value) -> double
	{
		return 1;
	};
	auto f2 = [](double &x_value) -> double
	{
		return 2.f + sqrt(x_value);
	};
	auto q1 = [](double &x_value) -> double
	{
		return 1;
	};
	auto q2 = [](double &x_value) -> double
	{
		return x_value * x_value;
	};


	std::vector<boundaryProblemPoint> points;
	points.push_back(p0);
	double step = p0.step;
	const double function_coef = 1.f / n;
	int i = 1;
	boundaryProblemPoint current_point;


	while (i < n)
	{
		current_point.number = i;
		current_point.step = step;
		current_point.x = points[i - 1].x + step;

		if ((current_point.x + step) < gap_point || current_point.x > gap_point)
		{
			current_point.a = function_coef * calc_a_value(points[i].x, points[i].x + step, k1);
			current_point.phi = function_coef * calc_phi_value(points[i].x - 0.5f, points[i].x + 0.5f, f1);
			current_point.d = function_coef * calc_phi_value(points[i].x - 0.5f, points[i].x + 0.5f, q1);
		}
		else if (1)
		{
			// добавить все случаи
		}

	return points;
}


inline std::vector<PPaint1> StartMethodBasic2(int n, PPaint1 p0, double control, double range, double xmax,  double a, double b)
{
	std::vector<PPaint1> points;
	
	points.push_back(p0);
	double h = p0.h;
	int i = 1;
	double xtmp;
	PPaint1 p;
	while (i < n + 1)
	{
		if (points[i - 1].xn > xmax - range)
			break;
		p.num = i;
		p.h = h;
		p.xn = points[i - 1].xn + h;
		if (p.xn > xmax)
		{
			h = h / 2;
			continue;
		}
		p.v1 = RK4_1(points[i - 1].xn, points[i - 1].v1, points[i - 1].v2, h);
		p.v2 = RK4_2(points[i - 1].xn, points[i - 1].v1, points[i - 1].v2, h, a, b);
		xtmp = p.xn - h / 2.0;
		p.v115 = RK4_1(points[i - 1].xn, points[i - 1].v1, points[i - 1].v2, h / 2.0);
		p.v215 = RK4_2(points[i - 1].xn, points[i - 1].v1, points[i - 1].v2, h / 2.0, a, b);
		p.v21 = RK4_1(xtmp, p.v115, p.v215, h / 2.0);
		p.v22 = RK4_2(xtmp, p.v115, p.v215, h / 2.0, a, b);
		p.S = (max(abs(p.v21), abs(p.v22)) - max(abs(p.v1), abs(p.v2))) / 15.0;
		if (abs(p.S) > control)
		{
			p.sub = -2;
			h = h / 2.0;
			continue;
		}
		if (abs(p.S) < control / 32.0)
		{
			h = 2.0 * h;
			p.sub = 2;
			points.push_back(p);
			p.sub = 0;
			i++;
			continue;
		}
		points.push_back(p);
		p.sub = 0;
		i++;
	}
	return points;
}

inline std::vector<PPaint1> Basic2(int n, PPaint1 p0, double control, double range, double xmax, double a, double b)
{
	std::vector<PPaint1> points;
	points.push_back(p0);
	double h = p0.h;
	int i = 1;
	PPaint1 p;
	while (i < n + 1)
	{
		if (points[i - 1].xn > xmax - range)
			break;
		p.num = i;
		p.h = h;
		p.xn = points[i - 1].xn + h;
		if (p.xn > xmax)
		{
			h = h / 2;
			continue;
		}
		p.v1 = RK4_1(points[i - 1].xn, points[i - 1].v1, points[i - 1].v2, h);
		p.v2 = RK4_2(points[i - 1].xn, points[i - 1].v1, points[i - 1].v2, h, a, b);
		p.S = (max(abs(p.v21), abs(p.v22)) - max(abs(p.v1), abs(p.v2))) / 15.0;
		points.push_back(p);
		i++;
	}
	return points;
}
