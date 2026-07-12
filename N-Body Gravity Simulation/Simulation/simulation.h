#pragma once

#include <vector>
#include <random>

using namespace std;

struct Vec2 { float x, y; };
struct QuadNode;
class QuadNodePool;

class Simulation{
    private:
        const Vec2 windowDimensions;
        const Vec2 corner;
        const long long N;
        const float timestep;
        const float MIN_SIZE = 0.0f;
        const float softening = 25.0f;
        const float theta_squared = 0.75;
        vector<Vec2> acc;
        vector<Vec2> new_pos;
        vector<Vec2> new_vel;
        //QuadNodePool NodePool;
    public:
        Simulation(const Vec2& Dimensions, const Vec2& corner, long long ParticleCount, float timestep);
        ~Simulation();
        vector<Vec2> generate_positions();
        vector<vector<Vec2>> generateGalaxy(Vec2 core1, Vec2 core2, float M1, float M2);
        void initialise(vector<vector<Vec2>>& Config, vector<float>& mass);
        vector<vector<Vec2>> NewState(vector<vector<Vec2>>& Config, vector<float>& mass, unsigned int& out);

    private:
        Vec2 com(Vec2& pos1, float m1, Vec2& pos2, float m2);
        unsigned int build_tree(Vec2& coord, QuadNode* node, QuadNodePool& nodes, float m);
        Vec2 GetForce(Vec2& coord, QuadNode* node);
        Vec2 getForce(Vec2& coord, QuadNode* node);
        void UpdatePosition(Vec2& coord, Vec2& Force, Vec2& Velocity, Vec2& target);
        void UpdateVelocity(Vec2& Force, Vec2& Velocity, Vec2& target);
};
