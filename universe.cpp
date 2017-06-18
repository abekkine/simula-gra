#include "universe.h"

#include <iostream>
#include <cmath>

#include "body.h"
#include "util.h"

namespace universe {

    int NumBodies;
    double SpreadRadius;
    int MassDistribution;

    void add_mass_distribution() {
        switch(MassDistribution) {
            case D_ROTATING_SPHERICAL:
                body::add_rotating_spherical(NumBodies, SpreadRadius);
                break;
            case D_RANDOM_CUBICAL:
                body::add_random_cubic(NumBodies, SpreadRadius);
                break;
            case D_RANDOM_SPHERICAL:
            case D_UNIFORM_CUBICAL:
            case D_UNIFORM_CUBICAL_W_NOISE:
                puts("Not implemented yet!");
                ::exit(1);
                break;
            default:
                break;
        }
    }

    void init() {
        add_mass_distribution();
    }

    void update() {

        #pragma omp parallel for
        for (int i=0; i<body::list.size(); i++) {
            body::list[i].pressure = 0.0;
            body::list[i].bound = 0;
            #pragma omp parallel for
            for (int j=0; j<body::list.size(); j++) {
                if (i!=j) {
                    body::apply_gravity(body::list[i], body::list[j]);
                }
            }

            if (body::list[i].bound > body::max_bound) {
                body::max_bound_ix = i;
                body::max_bound = body::list[i].bound;
                std::cout << "Bound(" << body::max_bound << ")" << std::endl;
                body::remaining_cam_steps = body::max_cam_steps;
            }
        }

        if (body::max_bound_ix != -1) {
            if (body::list[body::max_bound_ix].bound != body::max_bound) {
                body::max_bound_ix = -1;
                body::max_bound = 0;
            }
        }
    }

    void render() {

        display::world_mode();

        body::render();
        body::update();
    }
}

