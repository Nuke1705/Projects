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
        const long long N;
        const float timestep;
        const float MIN_SIZE = 5.0f;
    public:
        Simulation(const Vec2& Dimensions, long long ParticleCount, float timestep);
        ~Simulation();
        vector<Vec2> generate_positions();
        vector<vector<Vec2>> NewState(vector<vector<Vec2>>& Config, unsigned int& out);

    private:
        Vec2 com(Vec2& pos1, float m1, Vec2& pos2, float m2);
        unsigned int build_tree(Vec2& coord, QuadNode* node, QuadNodePool& nodes);
        Vec2 GetForce(Vec2& coord, QuadNode* node);
        Vec2 UpdatePosition(Vec2& coord, Vec2& Force, Vec2& Velocity);
        Vec2 UpdateVelocity(Vec2& Force, Vec2& Velocity, Vec2& coord);
};