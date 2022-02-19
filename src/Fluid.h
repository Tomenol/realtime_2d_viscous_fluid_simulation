#pragma once

#include <math.h>

#define borderTypeDensity					0x00
#define borderTypeVerticalVelocity			0x01
#define borderTypeHorizontVelocity			0x02

// Macros
#define swapMatrices(x, y){ float** tmp = x; x = y; y = tmp; }

class Fluid
{
private:
	// state at time t + dt
	float** _vx;
	float** _vy;

	float** _dye;
	float** _density;

	// state at time t
	float** _prevVx;
	float** _prevVy;

	float** _prevDye;
	float** _prevDensity;

	// Other
	int _size;
	float _timeStep;
	int _reynolds;

	int _gridSize;

public:
	Fluid();
	Fluid(int _s, int _gs, float _vx, float _vy, int _re, float _dens, float _dt);

	void addPrevDensity(int i, int j, float _dens);
	void addDensity(int i, int j, float _dens);

	float getDensity(int i, int j);
	float getVelocity(int x, int y, int _component);

	void setBorder(int _borderType, float** _x, int _velComponent);
	void setBorderCylinder(int _bt, float** _x, int r, int _vc);
	void setVelocity(int _x, int _y, float _vx, float _vy);
	void setPreviousVelocity(int _x, int _y, float _vx, float _vy);
	
	void diffuse(float** _x, float** _x0, float _amount, int _boundaryType, int _vc);
	void advect(float** _x, float** _x0, float** _vx0, float** _vy0, int _boundaryType, int _vc);
	void project(int N, float** u, float** v, float** p, float** div);

	void add_source(float** x, float** s, float dt);

	void simulate(float _dt);
	void vel_step(int N, float** u, float** v, float** u0, float** v0, float visc, float dt);
};




