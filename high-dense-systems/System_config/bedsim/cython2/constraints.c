#include <stdio.h>
#include <math.h>

#include "cconstraints.h"


int A(double angle, double major, double minor, double* ret) {
	double mat[2][2];

	double pref1 = 1/(major*major) - 1/(minor * minor);
	double pref2 = 1/(major*major) + 1/(minor * minor);
	mat[0][0] = 0.5 * (pref1 *  cos(2.0*angle) + pref2 );
	mat[0][1] = 0.5 * (pref1 *  sin(2.0*angle) );
	mat[1][0] = 0.5 * (pref1 *  sin(2.0*angle) );
	mat[1][1] = 0.5 * (pref1 * -cos(2.0*angle) + pref2 );
	ret = *mat;

	return 0;
}



double E(double t, double x, double y, double px, double py, double vx, double vy, double angle, double angvel, double major, double minor, double grad[]) {
    double cp, cm, a11, a12, a22, E, rx, ry;
    double a11p, a12p, a22p, Ep, Et;
    double Ex[2];

    rx = px + vx * t;
    ry = py + vy * t;
    angle += angvel * t;

    cm = (1. / (major * major) - 1. / (minor * minor)) / 2.;
    cp = (1. / (major * major) + 1. / (minor * minor)) / 2.;
    a11 =  cm * cos(2. * angle) + cp;
    a12 =  cm * sin(2. * angle);
    a22 = -cm * cos(2. * angle) + cp;

    E = a11 * (x - rx)*(x - rx) + 2. * a12 * (x - rx) * (y - ry) + a22 * (y - ry)*(y - ry) - 1.;

	if(grad) {
		a11p = -2. * cm * sin(2. * angle);
		a12p =  2. * cm * cos(2. * angle);
		a22p =  2. * cm * sin(2. * angle);

		Ep = a11p * (x - rx)*(x - rx) + 2. * a12p * (x - rx) * (y - ry) + a22p * (y - ry)*(y - ry);
		Ex[0] = 2. * (a11 * (x - rx) + a12 * (y - ry));
		Ex[1] = 2. * (a12 * (x - rx) + a22 * (y - ry));
		//Et = Ep * angvel - dot(Ex, v);
		Et = Ep * angvel - (Ex[0] * vx + Ex[1] * vy);

		grad[0] = Et;
		grad[1] = Ex[0];
		grad[2] = Ex[1];
	}

    return E;
}
