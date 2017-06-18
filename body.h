#ifndef BODY_H_
#define BODY_H_

#include <vector>
#include <cstdint>

namespace body {

    const int N = 3;

    extern double Seed;
    extern double Mass;
    extern double TimeStep;
    extern double SpeedRange;
    extern double GravityConstant;
    extern double CollisionRadius;
    extern double TerminalAcceleration;

    extern uint32_t max_bound;
    extern int32_t max_bound_ix;
    extern int32_t max_cam_steps;
    extern int32_t remaining_cam_steps;

    struct Body {
        double m;

        double& px;
        double& py;
        double& pz;

        double x[N];
        double vx[N];
        double ax[N];
        double pressure;
        uint32_t bound;

        Body(double& x0, double& y0, double& z0) :
            px(x0), py(y0), pz(z0) {}
    };

    extern std::vector< Body > list;

    extern float thetaX;
    extern float thetaY;
    extern float thetaZ;

    void add_mass_distribution();
    void add_rotating_spherical(int num, double radius);
    void add_random_spherical(int num, double radius);
    void add_uniform_cubic(int num, double length);
    void add_random_cubic(int num, double length);
    void apply_gravity(Body& from, Body& to);
    void render();
    void update();

    void update_camera_path();
    void drawCubeFrame();
}

#endif // BODY_H_

