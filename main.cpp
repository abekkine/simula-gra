#include <iostream>

#include "version.h"
#include "display.h"
#include "body.h"
#include "universe.h"

int main() {

    std::cout << "Version : " << VERSION_STR << std::endl;

    display::WorldSize = 200.0; // 500.0;
    display::RecordDuration = 300;
    display::VideoFile = "video.avi";
    display::ScreenWidth = 1920;
    display::ScreenHeight = 1080;
    display::withLookAt = false;

    body::ThetaXSpeed = 0.125;
    body::ThetaYSpeed = 0.25;
    body::ThetaZSpeed = 0.05;

    body::IgnitionLimit = 20;

    body::GravityConstant = 0.5;
    body::ExplosionConstant = 2.5;
    body::Seed = 2.99281; //1.22498;
    body::TimeStep = 0.01;
    body::Mass = 100.0;
    body::CollisionRadius = 0.8;
    body::TerminalAcceleration = 400.0;
    body::SpeedRange = 5.0;
    universe::MassDistribution = universe::DistributionType::D_ROTATING_SPHERICAL;
    //universe::MassDistribution = universe::DistributionType::D_RANDOM_SPHERICAL;
    //universe::MassDistribution = universe::DistributionType::D_UNIFORM_CUBICAL_W_NOISE;
    universe::NumBodies = 1500;
    universe::SpreadRadius = 200.0;

    display::init();
    universe::init();

    while (not display::quit()) {
        display::clear();
        universe::render();
        universe::update();
        display::update();
    }

    display::exit();

    return 0;
}

