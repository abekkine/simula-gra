#include "display.h"

#include <iostream>

#include <opencv2/opencv.hpp>

#include "util.h"
#include "body.h"

namespace display {

    int ScreenWidth;
    int ScreenHeight;
    double WorldSize;
    int RecordDuration;
    std::string VideoFile;

    int recordLength_;
    int screenWidth_;
    int screenHeight_;

    bool withLookAt = true;
    bool enableFullScreen_;
    bool enableBlend_;
    bool enableAA_;

    GLFWwindow* window_;

    void init() {
        puts("display::init()");

        screenWidth_ = ScreenWidth;
        screenHeight_ = ScreenHeight;
        recordLength_ = 25 * RecordDuration;

        enableFullScreen_ = false;
        enableBlend_ = true;
        enableAA_ = true;

        init_device();
    }

    void init_device() {
        puts("display::init_device()");

        glfwSetErrorCallback(errorCallback);
        if (not glfwInit()) {
            std::cerr << "GLFW initialization failed!" << std::endl;
            throw;
        }
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
        glfwWindowHint(GLFW_RESIZABLE, GLFW_FALSE);

        if (enableFullScreen_) {
            window_ = glfwCreateWindow(
                screenWidth_, screenHeight_,
                "F3",
                glfwGetPrimaryMonitor(),
                NULL);
        } else {
            window_ = glfwCreateWindow(
                screenWidth_, screenHeight_,
                "F3",
                NULL,
                NULL);
        }

        if (not window_) {
            glfwTerminate();
            std::cerr << "Unable to create GLFW window!" << std::endl;
            throw;
        }

        glfwSetKeyCallback(window_, keyCallback);
        glfwSetFramebufferSizeCallback(window_, resizeCallback);
        glfwSetCursorPosCallback(window_, cursorPositionCallback);
        glfwSetMouseButtonCallback(window_, mouseButtonCallback);

        glfwMakeContextCurrent(window_);
        glfwSwapInterval(1);

        if (enableBlend_) {
            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        }

        if (enableAA_) {
            glEnable(GL_POLYGON_SMOOTH);
            glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);

            glEnable(GL_LINE_SMOOTH);
            glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
        }
    }

    bool quit() {
        //puts("display::quit()");

        return glfwWindowShouldClose(window_);
    }

    void clear() {
        //puts("display::clear()");

        glClear(GL_COLOR_BUFFER_BIT);
    }

    void update() {
        //puts("display::update()");

        glfwSwapBuffers(window_);
        glfwPollEvents();

        save_buffer();
    }

    void exit() {
        puts("display::exit()");

        glfwTerminate();
    }

    void world_mode() {
        //puts("display::world_mode()");

        double wsize = WorldSize;
        double hsize = WorldSize * ScreenHeight / ScreenWidth;

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        if (withLookAt) {
            gluPerspective(60.0, wsize / hsize, 0.1, 1000);
        } else {
            glOrtho(-wsize, wsize, -hsize, hsize, -1000.0, 1000.0);
        }
        glMatrixMode(GL_MODELVIEW);
    }

    void ui_mode() {
        //puts("display::ui_mode()");

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glOrtho(0, screenWidth_, 0, screenHeight_, -1.0, 1.0);
        glMatrixMode(GL_MODELVIEW);
    }

    void errorCallback(int error, const char* description) {
        std::cout << "GLFX Error : " << description << std::endl;
    }

    void resizeCallback(GLFWwindow* w, int width, int height) {
        screenWidth_ = width;
        screenHeight_ = height;
        glViewport(0, 0, screenWidth_, screenHeight_);
    }

    void keyCallback(GLFWwindow* w, int key, int scancode, int action, int mods) {
        if (action == GLFW_RELEASE && key == GLFW_KEY_ESCAPE) {
            std::cout << "...ESC released..." << std::endl;
            glfwSetWindowShouldClose(w, GLFW_TRUE);
        } else if (action == GLFW_PRESS) {

            switch (key) {
                case GLFW_KEY_Z:
                    body::thetaX = 0.0; break;
                case GLFW_KEY_X:
                    body::thetaY = 0.0; break;
                case GLFW_KEY_C:
                    body::thetaZ = 0.0; break;
                default:
                    break;
            }
        }
    }

    static void mouseButtonCallback(GLFWwindow* w, int button, int action, int mods) {
        //std::cout << "...mouse button event..." << std::endl;
    }

    void cursorPositionCallback(GLFWwindow* w, double x, double y) {
        //std::cout << "cursor @ " << x << ", " << y << std::endl;
    }

    int frame = 0;
    cv::VideoWriter writer;
    void save_buffer() {

        if (frame == -1) {
            return;
        }
        if (frame == 0) {
            writer.open(
                VideoFile.c_str(),
                cv::VideoWriter::fourcc('X', '2', '6', '4'),
                25.0f,
                cv::Size(screenWidth_, screenHeight_),
                true);
            frame++;
        }
        else if (frame < recordLength_) {
            cv::Mat pixels( screenHeight_, screenWidth_, CV_8UC3 );
            glReadPixels(0, 0, screenWidth_, screenHeight_, GL_RGB, GL_UNSIGNED_BYTE, pixels.data);
            cv::Mat cv_pixels( screenHeight_, screenWidth_, CV_8UC3 );
            #pragma omp parallel for
            for (int y=0; y<screenHeight_; y++) {
                for (int x=0; x<screenWidth_; x++) {
                    cv_pixels.at<cv::Vec3b>(y, x)[2]
                        = pixels.at<cv::Vec3b>(screenHeight_-y-1, x)[0];
                    cv_pixels.at<cv::Vec3b>(y, x)[1]
                        = pixels.at<cv::Vec3b>(screenHeight_-y-1, x)[1];
                    cv_pixels.at<cv::Vec3b>(y, x)[0]
                        = pixels.at<cv::Vec3b>(screenHeight_-y-1, x)[2];
                }
            }
            writer << cv_pixels;
            frame++;
            if (frame % 100 == 0) {
                printf("%d/%d\n", frame, recordLength_);
            }
        }
        else if (frame == recordLength_) {
            writer.release();
            printf("Video is ready!\n");
            frame = -1;
            ::exit(0);
        }
    }
}
