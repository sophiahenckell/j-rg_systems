#include <stdio.h>
#include <math.h>


int project_on(double proj_angle, double angle, double major, double minor, double tnpoint[]) {
	/*! Project Ellipse with 'major' and 'minor' semi-axes on a line under 'proj_angle'. The ellipse angle is denoted as 'angle'.
	@param proj_angle Projection angle
	@param angle Angle of the ellipse
	@param major Semi major axis of the ellipse
	@param minor Semi minor axis of the ellipse
	@param tnpoint Empty array for output
	*/
	//double th = (proj_angle - angle) % (2 * M_PI);
	double th = fmod((proj_angle - angle),  (2 * M_PI));
	double direc_x = cos(th);
	double direc_y = sin(th);
	double xproj2, yproj2;
	double ppoint_x, ppoint_y, tnpoint_x, tnpoint_y, tnpref;

	if(tan(th) == 0.0) {
		xproj2 = major*major;
	} else {
		double m = -1.0/tan(th);
		if(m != m) { // m is NaN
			xproj2 = major*major;
		} else {
			xproj2 = (m * major) * (m * major) / (m*m + (minor*minor)/(major*major));
		}
	}
	yproj2 = minor*minor * (1.0 - xproj2/(major*major));
	ppoint_x = copysign(1.0, cos(th)) * sqrt(xproj2);
	ppoint_y = copysign(1.0, sin(th)) * sqrt(yproj2);

	tnpref = (ppoint_x * direc_x) + (ppoint_y * direc_y); // this is the dot(ppoint, direc) prefactor
	tnpoint_x = tnpref * direc_x; // /norm(direc)^2 not needed because direc has always length 1 here
	tnpoint_y = tnpref * direc_y;
	tnpoint[0] = tnpoint_x;
	tnpoint[1] = tnpoint_y;

	//rot = lambda th: np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
	
	//direc = np.dot(rot(th), np.array([1,0]))
	
	/*if np.tan(th) == 0: # catch zero division
		xproj2 = self.major**2
	else:
		m = -1/np.tan(th)
		if isfinite(m): # everything is just fine => this is the regular case
			xproj2 = (m * self.major)**2/(m**2 + (self.minor/self.major)**2)
		else: # catch infinite m's
			xproj2 = self.major**2
	
	yproj2 = self.minor**2 * (1-xproj2/self.major**2)
	*/
	
	//ppoint = np.sign(np.array([np.cos(th), np.sin(th)])) * np.sqrt([xproj2, yproj2])
	//tnpoint = np.dot(ppoint, direc) * direc/np.linalg.norm(direc)**2

	return 0;
}



int get_normal_projection_interval(double t, double proj_angle, double angle, double angvel, double major, double minor, double position[], double velocity[], double interval[]) {
	/*
	@param t project ellipse at time 't'
	@param proj_angle Angle of the _normal_ axis to project on
	*/
	double e_nsa, e_c;
	double nsa_dir_x = cos(proj_angle); 
	double nsa_dir_y = sin(proj_angle);
	double tnpoint[2];
	project_on(proj_angle - angvel * t, angle, major, minor, tnpoint);
	e_nsa = sqrt(tnpoint[0]*tnpoint[0] + tnpoint[1]*tnpoint[1]);
	e_c = (nsa_dir_x * position[0] + nsa_dir_y * position[1]) + (nsa_dir_x * velocity[0] + nsa_dir_y * velocity[1]) * t;
	interval[0] = e_c - e_nsa;
	interval[1] = e_c + e_nsa;

	//rot = lambda th: np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
	//nsa_dir = np.dot( rot(angle), np.array([1,0]) )
	
	//tnpoint = self._project_on(angle - self.angvel * t)
	//e_nsa = np.linalg.norm(tnpoint)
	//e_c = np.dot(nsa_dir, self.position) + np.dot(nsa_dir, self.velocity * t)
	//return (e_c - e_nsa, e_c + e_nsa)

	return 0;
}
