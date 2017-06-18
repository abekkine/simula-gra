#include "body.h"

#include <stdlib.h>
#include <math.h>

#include <iostream>
#include <random>

#include "display.h"
#include "util.h"

namespace body {

    double cx = 0.0;
    double cy = 0.0;
    double cz = 0.0;

    double Seed;
    double Mass;
    double TimeStep;
    double SpeedRange;
    double GravityConstant;
    double CollisionRadius;
    double TerminalAcceleration;

    std::vector< Body > list;
    int numBodies;
    double* positions;
    float* colors;

    uint32_t max_bound = 0;
    int32_t max_bound_ix = -1;
    int32_t max_cam_steps = 50;
    int32_t remaining_cam_steps = max_cam_steps;

    float cube_frame[] = {
         1,  1, -1,     -1,  1, -1,
        -1,  1, -1,     -1, -1, -1,
        -1, -1, -1,      1, -1, -1,
         1, -1, -1,      1,  1, -1,

         1,  1,  1,      1,  1, -1,
        -1,  1,  1,     -1,  1, -1,
         1, -1,  1,      1, -1, -1,
        -1, -1,  1,     -1, -1, -1,

         1,  1,  1,      1, -1,  1,
        -1,  1,  1,      1,  1,  1,
        -1, -1,  1,     -1,  1,  1,
         1, -1,  1,     -1, -1,  1
    };

    void add_rotating_spherical(int num, double radius) {

        numBodies = num;
        positions = new double[num * 3];
        colors = new float[num * 4];

        std::uniform_real_distribution<double> rdist(0.3 * radius, radius);
        std::uniform_real_distribution<double> pdist(-M_PI, M_PI);
        std::uniform_real_distribution<double> qdist(0.0, 2.0*M_PI);
        std::uniform_real_distribution<double> vdist(20.0, 40.0);
        std::uniform_real_distribution<double> wdist(-radius, radius);

        std::mt19937 rng;
        rng.seed(Seed);

        double r, p, q, v;
        for (int i=0; i<num; i++) {
            Body b(positions[3*i + 0], positions[3*i +1], positions[3*i + 2]);
            b.m = Mass;

            colors[4*i + 0] = 1.0;
            colors[4*i + 1] = 0.0;
            colors[4*i + 2] = 0.0;
            colors[4*i + 3] = 1.0;

            r = rdist(rng);
            p = pdist(rng);
            q = qdist(rng);
            b.x[0] = r * cos(p) * cos(q);
            b.x[1] = r * cos(p) * sin(q);
            b.x[2] = r * sin(p);

            b.px = b.x[0];
            b.py = b.x[1];
            b.pz = b.x[2];

            v = vdist(rng);
            b.vx[0] =  1.0 * b.x[1] * v / r;
            b.vx[1] = -1.0 * b.x[0] * v / r;
            b.vx[2] = 0.0;

            b.pressure = 0.0;
            b.bound = 0;

            list.push_back(b);
        }
    }

    void add_random_cubic(int num, double length) {

        numBodies = num;
        positions = new double[num * 3];

        std::uniform_real_distribution<double> xdist(-length, length);
        std::uniform_real_distribution<double> vdist(-SpeedRange, SpeedRange);

        std::mt19937 rng;
        rng.seed(std::random_device{}());
        for (int i=0; i<num; i++) {
            Body b(positions[3*i + 0], positions[3*i + 1], positions[3*i + 2]);
            b.m = Mass;

            for (int k=0; k<N; k++) {
                b.x[k] = xdist(rng);
                b.vx[k] = vdist(rng);
            }
            b.px = b.x[0];
            b.py = b.x[1];
            b.pz = b.x[2];

            b.pressure = 0.0;
            b.bound = 0;

            list.push_back(b);
        }
    }

    void add_random_spherical(int num, double radius) {
    }

    void add_uniform_cubic(int num, double length) {
    }

    float thetaX = 0.0;
    float thetaY = 0.0;
    float thetaZ = 0.0;
    void render() {

        glPushMatrix();
        glLoadIdentity();

        glRotatef(thetaX, 1.0, 0.0, 0.0);
        glRotatef(thetaY, 0.0, 1.0, 0.0);
        glRotatef(thetaZ, 0.0, 0.0, 1.0);

        // Focus on max bound index
        glTranslatef(-cx, -cy, -cz);
        //  End

        thetaX += 0.125;
        if (thetaX > 360.0) thetaX -= 360.0;
        thetaY += 0.25;
        if (thetaY > 360.0) thetaY -= 360.0;
        thetaZ += 0.05;
        if (thetaZ > 360.0) thetaZ -= 360.0;

        glEnableClientState(GL_VERTEX_ARRAY);
        glEnableClientState(GL_COLOR_ARRAY);
        glVertexPointer(3, GL_DOUBLE, 0, positions);
        glColorPointer(4, GL_FLOAT, 0, colors);

        glDrawArrays(GL_POINTS, 0, numBodies);
        glDisableClientState(GL_COLOR_ARRAY);
        glDisableClientState(GL_VERTEX_ARRAY);

        drawCubeFrame();

        // Draw a marker for given index
//        float mx = positions[2 * marker_ix + 0];
//        float my = positions[2 * marker_ix + 1];
//        glColor4f(1.0, 0.0, 0.0, 0.5);
//        glBegin(GL_LINE_LOOP);
//            glVertex2f(mx - 1.0, my - 1.0);
//            glVertex2f(mx - 1.0, my + 1.0);
//            glVertex2f(mx + 1.0, my + 1.0);
//            glVertex2f(mx + 1.0, my - 1.0);
//        glEnd();

        glPopMatrix();
    }

    void drawCubeFrame() {
       
        if (max_bound_ix == -1) return;

        glPushMatrix();

        glTranslatef(list[max_bound_ix].px,
            list[max_bound_ix].py,
            list[max_bound_ix].pz);

        glEnableClientState(GL_VERTEX_ARRAY);
        glVertexPointer(3, GL_FLOAT, 0, cube_frame);

        glScalef(2.0, 2.0, 2.0);
        glColor4f(1.0, 0.0, 0.0, 0.5);
        glDrawArrays(GL_LINES, 0, 24);

        glDisableClientState(GL_VERTEX_ARRAY);

        glPopMatrix();
    }

    void update() {

        double timeStep = TimeStep;

        //util::start_clock();
        #pragma omp parallel for
        for (int i=0; i<numBodies; i++) {

            for(int k=0; k<N; k++) {

                list[i].vx[k] += list[i].ax[k] * timeStep;
                list[i].x[k] += list[i].vx[k] * timeStep;
                list[i].ax[k] = 0.0;
            }

            list[i].px = list[i].x[0];
            list[i].py = list[i].x[1];
            list[i].pz = list[i].x[2];

            // Update colors by pressure.
            const float pmax = 4.5;
            const float pmin = -2.0;
            float c = std::log10(list[i].pressure);
            if (c < pmin) c = pmin;
            else if (c > pmax) c = pmax;
            c -= pmin;
            c = c / (pmax - pmin);
            colors[4*i + 0] = c;
            colors[4*i + 1] = c;
            colors[4*i + 2] = c;
            colors[4*i + 3] = 1;
        }
        //util::stop_clock();

        update_camera_path();
    }

    void update_camera_path() {

        if (max_bound_ix != -1) {
            float dx = list[max_bound_ix].px - cx;
            float dy = list[max_bound_ix].py - cy;
            float dz = list[max_bound_ix].pz - cz;

            if (remaining_cam_steps > 2) {
                cx += dx / remaining_cam_steps;
                cy += dy / remaining_cam_steps;
                cz += dz / remaining_cam_steps;
                remaining_cam_steps--;
            }
            else {
                cx += 0.1 * dx;
                cy += 0.1 * dy;
                cz += 0.1 * dz;
            }
            
        }
    }

    void apply_gravity(Body& from, Body& to) {
        const double K = GravityConstant;
        double dx[N];
        double dr, dr2;
        double to_ax[N], to_a;
        double Kmr;
        double vx[N];

        dr2 = 0.0;
        for (int k=0; k<N; k++) {
            dx[k] = to.x[k] - from.x[k];
            dr2 += dx[k] * dx[k];
        }
        dr = sqrt(dr2);

        to_a = 0.0;
        double drN = 1.0;
        for (int k=0; k<N; k++) drN *= dr;
        Kmr = K * to.m / drN;
        for (int k=0; k<N; k++) {
            to_ax[k] = Kmr * dx[k];
            to_a += to_ax[k] * to_ax[k];
        }

//        if (to_a > 400.0) {
        if (dr < CollisionRadius && to_a > TerminalAcceleration) {

            for (int k=0; k<N; k++) {
                vx[k] = (from.m * from.vx[k] + to.m * to.vx[k]) / (from.m + to.m);
                from.vx[k] = to.vx[k] = vx[k];
            }
        
            from.pressure += TerminalAcceleration;
            from.bound++;

        }
        else {
            for (int k=0; k<N; k++) {
                from.ax[k] += to_ax[k];
            }

            from.pressure += to_a;
        }
    }
}

