#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "res\headers\renderer.h"
#include "res\headers\vertexbuffer.h"
#include "res\headers\indexbuffer.h"
#include "res\headers\vertexarray.h"
#include "res\headers\shader.h"
#include "res\headers\VertexBufferLayout.h"
#include "res\headers\window.hpp"

#include "res\glm\glm.hpp"
#include "res\glm\gtc\matrix_transform.hpp"
#include "res\imgui\imgui.h"
#include "res\imgui\imgui_impl_glfw.h"
#include "res\imgui\imgui_impl_opengl3.h"

#include "Simulation\simulation.h"

#include <iostream>
#include <string>
#include <vector>
#include <array>

using namespace std;
using namespace glm;
using namespace Engine;
//-----------------------------------------------------------------------------------
#define WIDTH 1920.0f
#define HEIGHT 1080.0f
#define ParticleCount 20000
#define timestep 0.1f
//-----------------------------------------------------------------------------------

vector<float> ParticleObjectPoints(vector<Vec2>& positions){

    float sign[4] = {0.0f, 0.0f, 1.0f, 1.0f};
    vector<float> result(8*ParticleCount,0);

    for(int i = 0; i < ParticleCount; ++i){
        for(int j = 0; j < 4; ++j){
            result[2*j + 8*i] = positions[i].x + (2*sign[(j+1)%4]-1);
            result[2*j + 1 + 8*i] = positions[i].y + (2*sign[j]-1);
        }

    }
    return result;
}
// 0 1 2 3
// - + + -
// - - + +

vector<unsigned int> ParticleIndices(){

    unsigned int reference[6] = {
        0, 1, 2,
        2, 3, 0
    };
    vector<unsigned int> result(6*ParticleCount,0);

    for(int i = 0; i < ParticleCount; ++i){
        for(int j = 0; j < 6; ++j){
            result[6*i + j] = 4*i + reference[j];
        }
    }
    return result;
}
//-----------------------------------------------------------------------------------
bool paused = false;
int main(void){
    Window window(WIDTH , HEIGHT, "N-Body");
    Simulation gravity({WIDTH, HEIGHT}, ParticleCount, timestep);

    {
    
    vector<Vec2> original_position = gravity.generate_positions();
    vector<Vec2> Velocities(ParticleCount, {0,0});
    vector<vector<Vec2>> Config = {original_position, Velocities};
    vector<unsigned int> indices = ParticleIndices();
//-----------------------------------------------------------------------------------
    VertexBuffer vb(nullptr, 8*ParticleCount*sizeof(float), true);
    VertexBufferLayout layout;
    layout.Push(GL_FLOAT, 2);
    VertexArray va;
    va.AddBuffer(vb, layout);

    IndexBuffer ib(indices.data(), 6*ParticleCount);

    mat4 proj = ortho(0.0f, WIDTH, 0.0f, HEIGHT, -1.0f, 1.0f);
    mat4 view = translate(mat4(1.0f), vec3(0,0,0));
    mat4 model = translate(mat4(1.0f), vec3(0,0,0));

    mat4 mvp = proj * view * model;

    Shader shader("res/shaders/simple.shader");

    Renderer renderer;

    shader.Unbind();
    //vb.Unbind();
    //va.Unbind();
    ib.Unbind();

//-----------------------------------------------------------------------------------
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    //io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls
    //io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;      // Enable Gamepad Controls
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(window.getNativeWindow(), true);
    ImGui_ImplOpenGL3_Init("#version 130");
//-----------------------------------------------------------------------------------
    bool show_demo_window = true;
    bool show_another_window = false;
    ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);
//-----------------------------------------------------------------------------------
    unsigned int nodeCount = 0;
    /* Loop until the user closes the window */
    while (!window.shouldClose())
    {
        /* Render here */
        renderer.Clear();
//-----------------------------------------------------------------------------------
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
//-----------------------------------------------------------------------------------        
        //glDrawArrays(GL_TRIANGLES, 0, 6);
        vector<float> particles = ParticleObjectPoints(Config[0]);
        shader.Bind();
        //shader.SetUniform1f("u_time", renderer.GetTime());
        shader.SetUniformMat4f("u_MVP", mvp);
        
        vb.Update(particles.data(), particles.size()*sizeof(float));
        renderer.Draw(va, ib, shader);
        if(!paused){
            Config = gravity.NewState(Config, nodeCount);
        }   
//-----------------------------------------------------------------------------------
       {

            ImGui::Begin("Simulation Data");                          // Create a window called "Hello, world!" and append into it.

            //ImGui::Text("This is some useful text.");               // Display some text (you can use a format strings too)

            //ImGui::SameLine();
            ImGui::Text("Particle Count = %d", ParticleCount);
            if (ImGui::Button("Pause")){                            // Buttons return true when clicked (most widgets return true when edited/activated)
                paused = !paused;
            }
            ImGui::SameLine();
            if (ImGui::Button("Reset")){                            // Buttons return true when clicked (most widgets return true when edited/activated)
                Config = {original_position, Velocities};
            }
            ImGui::Text("Time step = %f", timestep);
            ImGui::Text("Number of nodes = %d", nodeCount);
            ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / io.Framerate, io.Framerate);
            ImGui::End();
        }
       ImGui::Render();
       ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        /* Swap front and back buffers */
        window.swapBuffers();

        /* Poll for and process events */
        window.pollEvents();
    }
    }
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    return 0;
}