#ifndef _SOLVER_H_
#define _SOLVER_H_

#define XY_TO_ARRAY(i,j) ((i)+(N+2)*(j))
#define FOR_EACH_CELL for ( i=1 ; i<=N ; i++ ) { for ( j=1 ; j<=N ; j++ ) {
#define END_FOR }}
#define SWAP(x0,x) {float * tmp=x0;x0=x;x=tmp;}

class Solver
{
	float dt, diff, visc;
	unsigned N;
	float * u_prev, *v_prev, *dens_prev;

public:
	float * u, * v, * dens;
	void Init(unsigned N, float dt, float diff, float visc);
	void FreeData(void);
	void ClearData(void);
	bool AllocateData(void);
	void ClearPrevData(void);
	void AddDensity(unsigned i, unsigned j, float source);
	void AddVelocity(unsigned i, unsigned j, float forceX, float forceY);
	void Solve(void);
private:
	void DensStep(void);
	void VelStep(void);

	void AddSource(float * x, float * s);
	void SetBounds(int b, float * x);
	void LinSolve( int b, float * x, float * x0, float a, float c);
	void Diffuse(int b, float * x, float * x0);
	void Advect(int b, float * d, float * d0, float * u, float * v);
	void Project(float * u, float * v, float * p, float * div);
};

#endif
