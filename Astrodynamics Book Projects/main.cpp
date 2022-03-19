#include <iostream>
#include <iomanip>
#include <vector>
using namespace std;

#include "glm/glm.hpp"
using namespace glm;

// Kepler project suggested in the "Fundamentals of Astrodynamics" book. Websites that were helpful are listed below
// https://ssd.jpl.nasa.gov/horizons/ --sample data recorded by nasa
// https://kyleniemeyer.github.io/space-systems-notes/orbital-mechanics/two-body-problems.html --specific to kepler problem
// http://www.nssc.ac.cn/wxzygx/weixin/201607/P020160718380095698873.pdf --specific to kepler problem

#define mu_earth 398600 // lower value than google due to km units
#define SERIES_COUNT 1000

float vectorMagnitude(vec3 v) {
	float sum = powf(v.x, 2.0f) + powf(v.y, 2.0f) + powf(v.z, 2.0f);
	return sqrtf(sum);
}
int factorial(int n) {
	int result = n;
	for (int i = n - 1; i > 0; i--) {
		result = result * i;
	}
	return result;
}

float C(float z) {
	float result = 0.0f;

	if (z > 0.0f) {
		result = (1.0f - cosf(sqrt(z))) / z;
	}
	if (z == 0.0f) {
		result = 1.0f / 2.0f;
	}
	if (z < 0.0f) {
		result = (coshf(sqrtf(-z)) - 1) / (-z);
	}

	return result;
}
float S(float z) {
	float result = 0.0f;

	if (z > 0.0f) {
		result = (sqrtf(z) - sin(sqrtf(z))) / powf(sqrt(z), 3.0f);
	}
	if (z == 0.0f) {
		result = 1.0f / 6.0f;
	}
	if (z < 0.0f) {
		result = (sinh(sqrtf(-z)) - sqrtf(-z)) / powf(sqrtf(-z), 3.0f);
	}

	return result;
}

float xToT(float x, float a, vec3 r, vec3 v, float dt) {
	// p197 equation
	float z = powf(x, 2.0f) / a;
	float rMag = vectorMagnitude(r);
	float vMag = dot(r, v) / rMag;

	float termOne = ((rMag * vMag) / sqrtf(mu_earth)) * powf(x, 2.0f) * C(z);
	float termTwo = (1.0f - (rMag / a)) * powf(x, 3.0f) * S(z);
	float termThree = rMag * x;
	float termFour = -sqrtf(mu_earth) * dt;
	
	return termOne + termTwo + termThree + termFour;
}
float xToSlope(float x, float a, vec3 r, vec3 v) {
	float z = powf(x, 2.0f) / a;
	float rMag = vectorMagnitude(r);
	float vMag = dot(r, v) / rMag;

	float termOne = ((rMag * vMag) / sqrtf(mu_earth)) * x * (1.0f - z * S(z));
	float termTwo = (1.0f - (rMag / a)) * powf(x, 2.0f) * C(z);
	float termThree = rMag;

	return termOne + termTwo + termThree;
}

float fFunction(float x, vec3 r, float a) {
	float z = powf(x, 2.0f) / a;
	float rMag = vectorMagnitude(r);

	return 1.0f - (powf(x, 2.0f) / rMag) * C(z);
}
float gFunction(float x, float dt, float a) {
	float z = powf(x, 2.0f) / a;

	return dt - (powf(x, 3.0f) / sqrtf(mu_earth)) * S(z);
}

float fDotFunction(vec3 r0, vec3 r, float x, float a) {
	float r0Mag = vectorMagnitude(r0);
	float rMag = vectorMagnitude(r);

	float z = powf(x, 2.0f) / a;

	float result = sqrtf(mu_earth) / (rMag * r0Mag);
	result = result * ((1.0f / a) * powf(x, 3.0f) * S(z) - x);
	
	return result;
}
float gDotFunction(float x, float a, vec3 r) {
	float rMag = vectorMagnitude(r);
	float z = powf(x, 2.0f) / a;

	return 1.0f - (powf(x, 2.0f) / rMag) * C(z);
}

void projectKepler(vec3 r, vec3 v, float dt) {
	// Step one - calulate rMag and semi major axis
	float rMag = vectorMagnitude(r);
	float vMag = vectorMagnitude(v);

	float numerator = 1.0f;
	float denominator = (2.0f / rMag) - (powf(vMag, 2.0f) / mu_earth);
	
	float a = numerator / denominator;

	// Side variable - eccentricity
	vec3 eccentricityVector = (1.0f / mu_earth) * ((powf(vMag, 2.0f) - (mu_earth / rMag)) * r - dot(r, v) * v);
	
	// calulate universal variable x
	float currentX = sqrtf(mu_earth) * abs(1.0f / a) * dt; // Chobotov approximation
	float lastError = 0.0f;

	for (int i = 0; i < 100; i++) {
		float newPrediction = xToT(currentX, a, r, v, dt);
		float slope = xToSlope(currentX, a, r, v);

		lastError = newPrediction / slope;

		currentX = currentX - (newPrediction / slope);
	}

	// find r 
	float f = fFunction(currentX, r, a);
	float g = gFunction(currentX, dt, a);
	
	vec3 resultantRadius = f * r + g * v;
	
	// find v
	float fDot = fDotFunction(r, resultantRadius, currentX, a);
	float gDot = gDotFunction(currentX, a, resultantRadius);

	vec3 resultantVelocity = fDot * r + gDot * v;
	
	// Outputs
	cout << "Orbital Eccentricity: " << vectorMagnitude(eccentricityVector) << endl;
	cout << "Orbital Semi-Major Axis: " << a << endl << endl;

	cout << "Final Error (Universal Variable): " << lastError << ", Universal Variable: " << currentX << endl << endl;

	cout << "Initial Radius Vector (km): (" << r.x << ", " << r.y << ", " << r.z << ")" << endl;
	cout << "Initial Velocity Vector (km/s): (" << v.x << ", " << v.y << ", " << v.z << ")" << endl;
	cout << "Time Difference Between Predictions (s): " << dt << endl << endl;

	cout << "Resultant Radius Vector: (" << resultantRadius.x << ", " << resultantRadius.y << ", " << resultantRadius.z << ")" << endl;
	cout << "Resultant Velocity Vector: (" << resultantVelocity.x << ", " << resultantVelocity.y << ", " << resultantVelocity.z << ")" << endl;
}

int main() {
	// Find data here
	// https://ssd.jpl.nasa.gov/horizons/

	vec3 r = vec3(-3.884740293102836E+03, -1.254121118685534E+03, 5.428547232601027E+03);
	vec3 v = vec3(1.907196707173440E+00, -7.412518691489030E+00, -3.372855620910822E-01);
	float dt = 86400.0f;

	projectKepler(r, v, dt);

	system("pause");
	return 0;
}