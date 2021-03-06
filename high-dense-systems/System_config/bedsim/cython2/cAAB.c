 
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "caAB.h" 
// generated by Mathematica 
/*Parameters 
aA - major semi-axis of ellipsoid A 
aB - major semi-axis of ellipsoid B 
bA - minor semi-axis of ellipsoid A 
bB - minor semi-axis of ellipsoid A 
cA - minor semi-axis of ellipsoid A 
cB - minor semi-axis of ellipsoid A 
rA0,rB0 - x coordinate of ellipsoid A,B 
rA1,rB1 - y coordinate of ellipsoid A,B 
rA2,rB2 - z coordinate of ellipsoid A,B 
qA0,qB0 - x coordinate of the quarternion for ellipsoid A,B 
qA1,qB1 - y coordinate of the quarternion for ellipsoid A,B 
qA2,qB2 - z coordinate of the quarternion for ellipsoid A,B 
qA3,qB3 - scalar part of the quarternion for ellipsoid A,B 
lambda0 - factor of the convex composite function 
*/ 
 
 
long double aAB0(long double aA, long double aB, long double bA, long double bB, long double cA, long double cB, long double qA0, long double qA1, long double qA2, long double qA3, long double qB0, long double qB1, long double qB2, long double qB3, long double rA0, long double rA1, long double rA2, long double rB0, long double rB1, long double rB2, long double lambda0) { 
     return (2.*(qB0*qB2 + qB1*qB3)*(0. + (2.*(qB0*qB2 + qB1*qB3))/cB) + 2.*(qB0*qB1 - qB2*qB3)*(0. + (2.*(qB0*qB1 - qB2*qB3))/bB) + 2.*(-0.5 + pow(qB0,2) + pow(qB3,2))*(0. + (2.*(-0.5 + pow(qB0,2) + pow(qB3,2)))/aB))*(-rA0 + rB0) + (2.*(qB1*qB2 - qB0*qB3)*(0. + (2.*(qB0*qB2 + qB1*qB3))/cB) + 2.*(-0.5 + pow(qB1,2) + pow(qB3,2))*(0. + (2.*(qB0*qB1 - qB2*qB3))/bB) + 2.*(qB0*qB1 + qB2*qB3)*(0. + (2.*(-0.5 + pow(qB0,2) + pow(qB3,2)))/aB))*(-rA1 + rB1) + (2.*(-0.5 + pow(qB2,2) + pow(qB3,2))*(0. + (2.*(qB0*qB2 + qB1*qB3))/cB) + 2.*(qB1*qB2 + qB0*qB3)*(0. + (2.*(qB0*qB1 - qB2*qB3))/bB) + 2.*(qB0*qB2 - qB1*qB3)*(0. + (2.*(-0.5 + pow(qB0,2) + pow(qB3,2)))/aB))*(-rA2 + rB2) ; 
} 
 
long double aAB1(long double aA, long double aB, long double bA, long double bB, long double cA, long double cB, long double qA0, long double qA1, long double qA2, long double qA3, long double qB0, long double qB1, long double qB2, long double qB3, long double rA0, long double rA1, long double rA2, long double rB0, long double rB1, long double rB2, long double lambda0) { 
     return (2.*(qB0*qB2 + qB1*qB3)*(0. + (2.*(qB1*qB2 - qB0*qB3))/cB) + 2.*(-0.5 + pow(qB0,2) + pow(qB3,2))*(0. + (2.*(qB0*qB1 + qB2*qB3))/aB) + 2.*(qB0*qB1 - qB2*qB3)*(0. + (2.*(-0.5 + pow(qB1,2) + pow(qB3,2)))/bB))*(-rA0 + rB0) + (2.*(qB1*qB2 - qB0*qB3)*(0. + (2.*(qB1*qB2 - qB0*qB3))/cB) + 2.*(qB0*qB1 + qB2*qB3)*(0. + (2.*(qB0*qB1 + qB2*qB3))/aB) + 2.*(-0.5 + pow(qB1,2) + pow(qB3,2))*(0. + (2.*(-0.5 + pow(qB1,2) + pow(qB3,2)))/bB))*(-rA1 + rB1) + (2.*(-0.5 + pow(qB2,2) + pow(qB3,2))*(0. + (2.*(qB1*qB2 - qB0*qB3))/cB) + 2.*(qB0*qB2 - qB1*qB3)*(0. + (2.*(qB0*qB1 + qB2*qB3))/aB) + 2.*(qB1*qB2 + qB0*qB3)*(0. + (2.*(-0.5 + pow(qB1,2) + pow(qB3,2)))/bB))*(-rA2 + rB2) ; 
} 
 
long double aAB2(long double aA, long double aB, long double bA, long double bB, long double cA, long double cB, long double qA0, long double qA1, long double qA2, long double qA3, long double qB0, long double qB1, long double qB2, long double qB3, long double rA0, long double rA1, long double rA2, long double rB0, long double rB1, long double rB2, long double lambda0) { 
     return (2.*(qB0*qB1 - qB2*qB3)*(0. + (2.*(qB1*qB2 + qB0*qB3))/bB) + 2.*(-0.5 + pow(qB0,2) + pow(qB3,2))*(0. + (2.*(qB0*qB2 - qB1*qB3))/aB) + 2.*(qB0*qB2 + qB1*qB3)*(0. + (2.*(-0.5 + pow(qB2,2) + pow(qB3,2)))/cB))*(-rA0 + rB0) + (2.*(-0.5 + pow(qB1,2) + pow(qB3,2))*(0. + (2.*(qB1*qB2 + qB0*qB3))/bB) + 2.*(qB0*qB1 + qB2*qB3)*(0. + (2.*(qB0*qB2 - qB1*qB3))/aB) + 2.*(qB1*qB2 - qB0*qB3)*(0. + (2.*(-0.5 + pow(qB2,2) + pow(qB3,2)))/cB))*(-rA1 + rB1) + (2.*(qB1*qB2 + qB0*qB3)*(0. + (2.*(qB1*qB2 + qB0*qB3))/bB) + 2.*(qB0*qB2 - qB1*qB3)*(0. + (2.*(qB0*qB2 - qB1*qB3))/aB) + 2.*(-0.5 + pow(qB2,2) + pow(qB3,2))*(0. + (2.*(-0.5 + pow(qB2,2) + pow(qB3,2)))/cB))*(-rA2 + rB2) ; 
} 
 
 
