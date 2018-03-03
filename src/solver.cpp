#include "solver.h"
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <cstring>
#include <cmath>

void Solver::Init(unsigned N, float dt, float diff, float visc)
{
	this->dt = dt;
	this->diff = diff;
	this->visc = visc;
	this->N = N;
}
/*
----------------------------------------------------------------------
free/clear/allocate simulation data
----------------------------------------------------------------------
*/
void Solver::FreeData(void)
{
	if (this->u != NULL) free(this->u);
	if (this->v != NULL) free(this->v);
	if (this->dens != NULL) free(this->dens);
	
	if (this->v_prev != NULL) free(this->v_prev);
	if (this->u_prev != NULL) free(this->u_prev);
	if (this->dens_prev != NULL) free(this->dens_prev);
}

void Solver::ClearData(void)
{
	const int sizeBytes = (this->N + 2) * (this->N + 2) * sizeof(float);

	memset(this->u, 0, sizeBytes);
	memset(this->v, 0, sizeBytes);
	memset(this->dens, 0, sizeBytes);

	memset(this->u_prev, 0, sizeBytes);
	memset(this->v_prev, 0, sizeBytes);
	memset(this->dens_prev, 0, sizeBytes);
}

bool Solver::AllocateData(void)
{
	const int sizeBytes = (this->N + 2) * (this->N + 2) * sizeof(float);

	this->v = (float *)malloc(sizeBytes);
	this->u = (float *)malloc(sizeBytes);
	this->dens = (float *)malloc(sizeBytes);

	this->v_prev = (float *)malloc(sizeBytes);
	this->u_prev = (float *)malloc(sizeBytes);
	this->dens_prev = (float *)malloc(sizeBytes);

	if (this->v_prev == NULL || this->u_prev == NULL || this->dens_prev == NULL
		|| this->v == NULL || this->u == NULL || this->dens == NULL) {
		return false;
	}

	this->ClearData();
	return true;
}

void Solver::ClearPrevData() 
{
	const int sizeBytes = (this->N + 2) * (this->N + 2) * sizeof(float);

	memset(this->u_prev, 0, sizeBytes);
	memset(this->v_prev, 0, sizeBytes);
	memset(this->dens_prev, 0, sizeBytes);
}

void Solver::AddDensity(unsigned x, unsigned y, float source)
{
	this->dens_prev[XY_TO_ARRAY(x, y)] = source;
}

void Solver::AddVelocity(unsigned x, unsigned y, float forceX, float forceY)
{
	this->u_prev[XY_TO_ARRAY(x, y)] = forceX;
	this->v_prev[XY_TO_ARRAY(x, y)] = forceY;
}

void Solver::Solve()
{
	VelStep();
	DensStep();
}

void Solver::DensStep()
{
	AddSource(dens, dens_prev);			//Adding input density (dens_prev) to final density (dens).
	SWAP(dens_prev, dens)				//Swapping matrixes, because we want save the next result in dens, not in dens_prev.
	Diffuse(0, dens, dens_prev);		//Writing result in dens because we made the swap before. bi = dens_prev. The initial trash in dens matrix, doesnt matter, because it converges anyways.
	SWAP(dens_prev, dens)				//Swapping matrixes, because we want save the next result in dens, not in dens_prev.
	Advect(0, dens, dens_prev, u, v);	//Advect phase, result in dens.
}

void Solver::VelStep()
{
	AddSource(u, u_prev);
	AddSource(v, v_prev);
	SWAP (u_prev,u)			
	SWAP (v_prev, v)
	Diffuse(1, u, u_prev);  
	Diffuse(2, v, v_prev); 
	Project(u, v, u_prev, v_prev);		//Mass conserving.
	SWAP (u_prev,u)			
	SWAP (v_prev,v)
	Advect(1, u, u_prev, u_prev, v_prev);
	Advect(2, v, v_prev, u_prev, v_prev);
	Project(u, v, u_prev, v_prev);		//Mass conserving.
}

void Solver::AddSource(float * base, float * source)
{
//TODO: Teniendo en cuenta dt (Delta Time), incrementar el array base con nuestro source. Esto sirve tanto para añadir las nuevas densidades como las nuevas fuerzas.
	int i, j;
	FOR_EACH_CELL
		base[XY_TO_ARRAY(i, j)] += source[XY_TO_ARRAY(i, j)] * this->dt;
	END_FOR
}


void Solver::SetBounds(int b, float * x)
{
/*TODO:
Input b: 0, 1 or 2.
	0: borders = same value than the inner value.
	1: x axis borders inverted, y axis equal.
	2: y axis borders inverted, x axis equal.
	Corner values allways are mean value between associated edges.
*/
	for (int i = 1; i < this->N + 1; ++i) {
		// fila arriba
		x[XY_TO_ARRAY(0, i)] = b == 1 ? -x[XY_TO_ARRAY(1, i)] : x[XY_TO_ARRAY(1, i)];

		// fila abajo
		x[XY_TO_ARRAY(N + 1, i)] = b == 1 ? -x[XY_TO_ARRAY(N, i)] : x[XY_TO_ARRAY(N, i)];

		// columna izq
		x[XY_TO_ARRAY(i, 0)] = b == 2 ? -x[XY_TO_ARRAY(i, 1)] : x[XY_TO_ARRAY(i, 1)];

		// columna dcha
		x[XY_TO_ARRAY(i, N + 1)] = b == 2 ? -x[XY_TO_ARRAY(i, N)] : x[XY_TO_ARRAY(i, N)];
    }
}

/*
https://www.youtube.com/watch?v=62_RUX_hrT4
https://es.wikipedia.org/wiki/M%C3%A9todo_de_Gauss-Seidel <- Solución de valores independientes.
Despreciando posibles valores de x no contiguos, se simplifica mucho. Mirar diapositivas y la solución de Gauss Seidel de términos independientes.
Gauss Seidel -> Matrix x and x0
*/
void Solver::LinSolve(int b, float * x, float * x0, float aij, float aii)
{
// Se recomienda usar FOR_EACH_CELL, END_FOR y XY_TO_ARRAY.
	int i, j;
	for (int k = 0; k < 25; ++k) {
		FOR_EACH_CELL
			x[XY_TO_ARRAY(i, j)] = (x0[XY_TO_ARRAY(i, j)] - aij * (x[XY_TO_ARRAY(i + 1, j)] + x[XY_TO_ARRAY(i - 1, j)] +
				x[XY_TO_ARRAY(i, j + 1)] + x[XY_TO_ARRAY(i, j - 1)])) / aii;
		END_FOR

		SetBounds(b, x);
	}
}

/*
Nuestra función de difusión solo debe resolver el sistema de ecuaciones simplificado a las celdas contiguas de la casilla que queremos resolver,
por lo que solo con la entrada de dos valores, debemos poder obtener el resultado.
*/
void Solver::Diffuse(int b, float * x, float * x0)
{
// Solo necesitaremos pasar dos parámetros a nuestro resolutor de sistemas de ecuaciones de Gauss Seidel. Calculamos dichos valores y llamamos a la resolución del sistema.
	float a = dt * diff * N * N;
	float aii = 1 + 4*a;
	float aij = -a;
	LinSolve(b, x, x0, aij, aii);
}

/*
d is overwrited with the initial d0 data and affected by the u & v vectorfield.
Hay que tener en cuenta que el centro de las casillas representa la posición entera dentro de la casilla, por lo que los bordes estan
en las posiciones x,5.
*/
void Solver::Advect(int b, float * d, float * d0, float * u, float * v)
{
	int i, j;
	FOR_EACH_CELL
		float pu = static_cast<float>(i) - u[XY_TO_ARRAY(i, j)] * dt * N;
	    float pv = static_cast<float>(j) - v[XY_TO_ARRAY(i, j)] * dt * N;

		pu = pu < 0.5f ? 0.5f : pu > N + 0.5f ? N + 0.5f : pu;
		pv = pv < 0.5f ? 0.5f : pv > N + 0.5f ? N + 0.5f : pv;
    
		// Para acceder a d0
	    int iu = static_cast<int>(pu);
	    int iv = static_cast<int>(pv);
    
		// diffs
	    float du = std::abs(pu - iu);
	    float dv = std::abs(pv - iv);
    
	    d[XY_TO_ARRAY(i, j)] = ((1 - du) * (1 - dv)) * d0[XY_TO_ARRAY(iu, iv)]
		                     + ((1 - du) * dv)       * d0[XY_TO_ARRAY(iu, iv + 1)]
		                     + (du       * (1 - dv)) * d0[XY_TO_ARRAY(iu + 1, iv)]
		                     + (du       * dv)       * d0[XY_TO_ARRAY(iu + 1, iv + 1)];
	END_FOR
}

/*
Se encarga de estabilizar el fluido y hacerlo conservativo de masa. Se usa solo en las matrices de velocidades.
No necesaria implementación por su complejidad.
*/
void Solver::Project(float * u, float * v, float * p, float * div)
{
	int i, j;

	FOR_EACH_CELL
		div[XY_TO_ARRAY(i, j)] = -0.5f*(u[XY_TO_ARRAY(i + 1, j)] - u[XY_TO_ARRAY(i - 1, j)] + v[XY_TO_ARRAY(i, j + 1)] - v[XY_TO_ARRAY(i, j - 1)]) / N;
		p[XY_TO_ARRAY(i, j)] = 0;
	END_FOR
	SetBounds(0, div);
	SetBounds(0, p);

	LinSolve(0, p, div, 1, 4);

	//Aproximamos: Laplaciano de q a su gradiente.
	FOR_EACH_CELL
		u[XY_TO_ARRAY(i, j)] -= 0.5f*N*(p[XY_TO_ARRAY(i + 1, j)] - p[XY_TO_ARRAY(i - 1, j)]);
		v[XY_TO_ARRAY(i, j)] -= 0.5f*N*(p[XY_TO_ARRAY(i, j + 1)] - p[XY_TO_ARRAY(i, j - 1)]);
	END_FOR
	SetBounds(1, u);
	SetBounds(2, v);
}