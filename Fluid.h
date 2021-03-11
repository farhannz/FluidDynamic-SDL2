#include <stdlib.h>
#include <math.h> 
#define IX(x, y) (constrain(x,0,N-1) + (constrain(y,0,N-1)*N))
class Fluid{
    private:
        int size;
        float dt;
        float diff;
        float visc;
        float *s;
        float *density;
        
        float *Vx;
        float *Vy;

        float *Vx0;
        float *Vy0;
    public:
        Fluid(int size, float diffusion, float viscosity, float dt);
        ~Fluid();
        int Size();
        float Density(int x,int y);
        void addDensity(int x, int y,float amount);
        void addVelocity(int x, int y, float amountX, float amountY);
        void step();
        void renderGrid();
};

void set_bnd(int bound, float *x , int N);
void lin_solve(int b, float *x, float *x0, float a, float c, int iter, int N);
void diffuse (int b, float* x, float* x0, float diff, float dt, int iter, int N);
void project(float *velocX, float *velocY, float *p, float *div, int iter, int N);
void advect(int b, float* d, float* d0, float* velocX, float* velocY, float dt, int N);
int constrain(int val, int low, int high);