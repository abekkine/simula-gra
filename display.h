#ifndef DISPLAY_H_
#define DISPLAY_H_

#include <GLFW/glfw3.h>

#include <string>

namespace display {

    extern int ScreenWidth;
    extern int ScreenHeight;
    extern double WorldSize;
    extern int RecordDuration;
    extern std::string VideoFile;

    void init();
    void exit();
    bool quit();
    void update();
    void clear();
    void init_device();
    void world_mode();
    void save_buffer();
    void ui_mode();

    static void errorCallback(int error, const char* description);
    static void resizeCallback(GLFWwindow* w, int width, int height);
    static void keyCallback(GLFWwindow* w, int key, int scancode, int action, int mods);
    static void mouseButtonCallback(GLFWwindow* w, int button, int action, int mods);
    static void cursorPositionCallback(GLFWwindow* w, double x, double y);
}

#endif // DISPLAY_H_

