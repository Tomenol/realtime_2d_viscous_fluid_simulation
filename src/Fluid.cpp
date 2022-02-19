#include "Fluid.h"

Fluid::Fluid()
{

}

Fluid::Fluid(int _s, int _gs, float _vx, float _vy, int _re, float _dens, float _dt) : _size(_s), _gridSize(_gs)
{
	this->_dye = new float* [this->_size];
	this->_density = new float* [this->_size];
	this->_vx = new float* [this->_size];
	this->_vy = new float* [this->_size];

	this->_prevDye = new float* [this->_size];
	this->_prevDensity = new float* [this->_size];
	this->_prevVx = new float* [this->_size];
	this->_prevVy = new float* [this->_size];

	this->_reynolds = _re;
	this->_timeStep = _dt;

	for (int i = 0; i < this->_size; i++)
	{
		_density[i] = new float[this->_size];
		this->_vx[i] = new float[this->_size];
		this->_vy[i] = new float[this->_size];

		this->_prevDensity[i] = new float[this->_size];

		this->_prevVx[i] = new float[this->_size];
		this->_prevVy[i] = new float[this->_size];

		for (int j = 0; j < _size; j++)
		{
			this->_density[i][j] = _dens;
			this->_prevDensity[i][j] = _dens;

			setVelocity(i, j, _vx, _vy);
			setPreviousVelocity(i, j, _vx, _vy);
		}
	}
}

float Fluid::getDensity(int i, int j)
{
	return this->_density[i][j];
}

void Fluid::addPrevDensity(int i, int j, float _dens)
{
	this->_prevDensity[i][j] = _dens;
}

void Fluid::addDensity(int i, int j, float _dens)
{
	this->_density[i][j] = _dens;
}

float Fluid::getVelocity(int x, int y, int _component)
{
	if (_component == 0) return this->_vx[x][y];
	if (_component == 1) return this->_vy[x][y];

	return -1;
}

void Fluid::setVelocity(int _x, int _y, float _vx, float _vy)
{
	this->_vx[_x][_y] = _vx;
	this->_vy[_x][_y] = _vy;
}

void Fluid::setPreviousVelocity(int _x, int _y, float _vx, float _vy)
{
	this->_prevVx[_x][_y] = _vx;
	this->_prevVy[_x][_y] = _vy;
}

void Fluid::diffuse(float** _x, float** _x0, float _amount, int _boundaryType, int _vc)
{
	float _diff = _amount * this->_timeStep * this->_size * this->_size;

	for (int k = 0; k < 20; k++)
	{
		for (int i = 1; i < this->_size - 2; i++)
			for (int j = 1; j < this->_size - 2; j++)
				_x[i][j] = (_x0[i][j] + _diff * ( _x[i - 1][j] + _x[i + 1][j] + _x[i][j - 1] + _x[i][j + 1])) / (1 + 4 * _diff);

		setBorder(_boundaryType, _x, _vc);
	}
}

void Fluid::setBorder(int _borderType, float** _x, int _velComponent)
{
	for (int i = 0; i < this->_size - 1; i++)
	{
		_x[0][i] = (_borderType == 1 ? -_x[1][i] : _x[1][i]);
		_x[this->_size - 1][i] = (_borderType == 1 ? -_x[this->_size - 2][i] : _x[this->_size - 2][i]);
		_x[i][0] = (_borderType == 2 ? -_x[i][1] : _x[i][1]);
		_x[i][this->_size - 1] = (_borderType == 2 ? -_x[i][this->_size - 2] : _x[i][this->_size - 2]);
	}

	_x[0][0] = 0.5 * ((float)_x[0][1] + (float)_x[1][0]);
	_x[this->_size - 1][0] = 0.5 * ((float)_x[this->_size - 2][0] + (float)_x[_size - 1][1]);
	_x[0][this->_size - 1] = 0.5 * ((float)_x[0][this->_size - 2] + (float)_x[1][_size - 1]);
	_x[this->_size - 1][_size - 1] = 0.5 * ((float)_x[this->_size - 2][this->_size - 1] + (float)_x[this->_size - 1][this->_size - 2]);

	//setBorderCylinder(_borderType, _x, 20, 0);
}

void Fluid::advect(float** _x, float** _x0, float** _vx0, float** _vy0, int _boundaryType, int _vc)
{
	for (int i = 1; i < this->_size - 1; i++)
		for (int j = 1; j < this->_size - 1; j++)
		{
			// get the previous position of the fluid particle (semi Lagrangian)
			float _prevX = i - _vx0[i][j] * this->_timeStep;
			float _prevY = j - _vy0[i][j] * this->_timeStep;

			// boundary conditions
			if (_prevX > this->_size - 2) _prevX = this->_size - 2;
			if (_prevX < 0) _prevX = 0;

			if (_prevY > this->_size - 2) _prevY = this->_size - 2;
			if (_prevY < 0) _prevY = 0;

			// get closest point on the grid
			int i0 = (int)(_prevX);
			int j0 = (int)(_prevY);

			float a1 = _prevX - i0;
			float a2 = _prevY - j0;

			// get previous values by interpolation
			_x[i][j] = (_x0[i0 + 1][j0] * (1 - a2) + _x0[i0 + 1][j0 + 1] * a2) * a1 + (_x0[i0][j0] * (1 - a2) + _x0[i0][j0 + 1] * a2) * (1 - a1);
		}
	setBorder(_boundaryType, _x, _vc);
}

void Fluid::project(int N, float** u, float** v, float** p, float** div)
{
	int i, j, k;
	float h;
	h = 1.0 / N;

	for (i = 1; i < N - 1; i++) {
		for (j = 1; j < N - 1; j++) {
			div[i][j] = 0.5 * h * (u[i + 1][j] - u[i - 1][j] + v[i][j + 1] - v[i][j - 1]); // calculate the divergence of the flow
			p[i][j] = 0; // set pressures to 0
		}
	}

	setBorder(0, div, 0);
	setBorder(0, p, 0);

	for (k = 0; k < 20; k++)
	{
		for (i = 1; i < N - 1; i++)
			for (j = 1; j < N - 1; j++)
				p[i][j] = (-div[i][j] + p[i - 1][j] + p[i + 1][j] + p[i][j - 1] + p[i][j + 1]) / 4;

		setBorder(0, p, 0);
	}

	// set speed vectors
	for (i = 1; i < N - 1; i++)
	{
		for (j = 1; j < N - 1; j++)
		{
			u[i][j] -= 0.5 * (p[i + 1][j] - p[i - 1][j]) / h;
			v[i][j] -= 0.5 * (p[i][j + 1] - p[i][j - 1]) / h;
		}
	}

	setBorder(1, u, 1);
	setBorder(2, v, 2);
}

void Fluid::add_source(float** x, float** s, float dt)
{
	int i, j;

	for (i = 0; i < this->_size; i++)
		for (j = 0; j < this->_size; j++)
		{
			x[i][j] += dt * s[i][j];
			//s[i][j] -= dt * s[i][j];
		}
}

void Fluid::simulate(float _dt)
{
	_timeStep = _dt;

	vel_step(_size, _vx, _vy, _prevVx, _prevVy, 1e-8, _timeStep);

	diffuse(this->_prevDensity, this->_density, 1e-4, 0, 0);
	advect(this->_density, this->_prevDensity, this->_vx, this->_vy, 0, 0);
}

void Fluid::vel_step(int N, float** u, float** v, float** u0, float** v0, float visc, float dt)
{
	add_source(u, u0, dt);
	add_source(v, v0, dt);

	diffuse(u0, u, visc, 1, 1);
	diffuse(v0, v, visc, 2, 2);

	project(N, u0, v0, u, v);

	advect(u, u0, u0, v0, 1, 1);
	advect(v, v0, u0, v0, 2, 2);

	project(N, u, v, u0, v0);
}
