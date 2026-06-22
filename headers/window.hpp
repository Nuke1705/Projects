#pragma once

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <string>

namespace Engine{

    class Window{
        private:
            GLFWwindow* m_window;
        public:
            Window(int width, int height, const std::string& name);

            ~Window();
            
            bool shouldClose() const;

            void pollEvents();

            void swapBuffers();

        GLFWwindow* getNativeWindow();
    };
}