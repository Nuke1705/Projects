//#include <SFML/Graphics.hpp>
#include <vector>
#include <random>
#include "simulation.h"
#include <omp.h>
#include <chrono>
#include <iostream>
//using namespace sf;
using namespace std;

//struct Vec2 { float x, y; };

struct QuadNode {
    // bounds of this node
    float x, y;        // center or top-left, depending on your convention
    float size_x;      // width of this square
    float size_y;
    // children
    QuadNode* children[4] = {nullptr, nullptr, nullptr, nullptr};

    // objects in this node
    Vec2 COM;
    float mass = 0;
    bool first = true;;
    vector<Vec2> same_node;
    QuadNode* nxtfree = nullptr;
    //node functions
    bool isleaf(){
        return (children[0]==nullptr)&&(children[1]==nullptr)&&(children[2]==nullptr)&&(children[3]==nullptr);
    }
};

class QuadNodePool{
public:
    QuadNodePool(size_t blockSize)
    : blockSize(blockSize) {}

    QuadNode* allocate(){
        if(!freeList) addBlock();

        QuadNode* n = freeList;

        freeList = freeList->nxtfree;

        fill(n->children, n->children + 4, nullptr);
        n->mass = 0;
        n->nxtfree = nullptr;
        n->COM = {0.0f, 0.0f};

        return n;
    }

    void deallocate(QuadNode* n) {
        n->nxtfree = freeList;
        freeList = n;
    }

    ~QuadNodePool() {
        for (auto block : blocks)
            delete[] block;
    }

private:
    size_t blockSize;
    vector<QuadNode*> blocks;
    QuadNode* freeList = nullptr;
    void addBlock() {
        QuadNode* block = new QuadNode[blockSize];
        blocks.push_back(block);

        // chain into freelist
        for (std::size_t i = 0; i < blockSize - 1; ++i)
            block[i].nxtfree = &block[i + 1];

        block[blockSize - 1].nxtfree = freeList;
        freeList = block;
    }

};

Simulation::Simulation(const Vec2& Dimensions, const Vec2& corner, long long ParticleCount, float timestep)
: windowDimensions(Dimensions), N(ParticleCount), timestep(timestep), corner(corner){
    acc.resize(N);
    new_pos.resize(N);
    new_vel.resize(N);
}

Simulation::~Simulation(){
}

Vec2 Simulation::com(Vec2& pos1, float m1, Vec2& pos2, float m2){
    Vec2 temp;
    float M = m1 + m2;
    temp.x = ((pos1.x*m1)+(pos2.x*m2))/M;
    temp.y = ((pos1.y*m1)+(pos2.y*m2))/M;
    return temp;
}

unsigned int Simulation::build_tree(Vec2& coord, QuadNode* node, QuadNodePool& nodes, float m){
    unsigned int nodecount = 0;
    if(node->mass == 0){
        node->mass = m;   //root case
        node->COM = coord;
        nodecount++;
        return nodecount;
    }
    if(node->isleaf()){ //if the node has a particle but no children
        if(min(node->size_x/2, node->size_y/2) < MIN_SIZE){
            if(node->first){
                node->same_node.push_back(node->COM);
                node->first = false;
            }
            node->same_node.push_back(coord);
            Vec2 ncom = com(coord, m, node->COM, node->mass);
            node->mass += m;
            node->COM = ncom; 
            return nodecount;
        }
        QuadNode* branch1 = nodes.allocate();
        branch1->COM = node->COM;
        int a = (node->COM.x >= node->x + node->size_x / 2);
        int b = (node->COM.y >= node->y + node->size_y / 2);
        int c = 2*b + a;
        node->children[c] = branch1;
        branch1->mass = node->mass;
        branch1->x = node->x + a*(node->size_x/2);
        branch1->y = node->y + b*(node->size_y/2);
        branch1->size_x = node->size_x/2;
        branch1->size_y = node->size_y/2;

        int A = (coord.x >= node->x + node->size_x / 2);
        int B = (coord.y >= node->y + node->size_y / 2);
        int C = 2*B + A;
            
        if(C == c){
            nodecount += build_tree(coord, branch1, nodes, m);
        }
        else{
            QuadNode* branch2 = nodes.allocate();
            branch2->COM = coord;
            node->children[C] = branch2;
            branch2->mass = m;
            branch2->x = node->x + A*(node->size_x/2);
            branch2->y = node->y + B*(node->size_y/2);
            branch2->size_x = node->size_x/2;
            branch2->size_y = node->size_y/2;
            nodecount++;
        }   
        Vec2 ncom = com(coord, m, node->COM, node->mass);
        node->COM = ncom;
        node->mass += m;
        nodecount++;
        return nodecount;
    }
    else{ //has atleast one child i.e. mass>1
        int A = (coord.x >= node->x + node->size_x / 2);
        int B = (coord.y >= node->y + node->size_y / 2);
        int C = 2*B + A;
        Vec2 ncom = com(coord, m, node->COM, node->mass);
        node->COM = ncom;
        node->mass += m;
        if(node->children[C] == nullptr){
            QuadNode* branch2 = nodes.allocate();
            branch2->COM = coord;
            node->children[C] = branch2;
            branch2->mass = m;
            branch2->x = node->x + A*(node->size_x/2);
            branch2->y = node->y + B*(node->size_y/2);
            branch2->size_x = node->size_x/2;
            branch2->size_y = node->size_y/2;
            nodecount++;
        }
        else{
            nodecount += build_tree(coord, node->children[C], nodes, m);
        }
        return nodecount;
    }
}

Vec2 Simulation::GetForce(Vec2& coord, QuadNode* node){ //for aggressive approximation (Not used)
    Vec2 force = {0.0f, 0.0f};
    if(node->isleaf()){
        for(const Vec2& pos : node->same_node){
            float x = pos.x;
            float y = pos.y;

            float dx = x - coord.x;
            float dy = y - coord.y;
            if(dx*dx + dy*dy > 1e-6f){
                float r2 = dx*dx + dy*dy + softening;//1e-4f; // softening
                float r = sqrt(r2);
                float inv_r3 = 1.0f / (r2 * r);

                float Fx = (dx) * inv_r3;
                float Fy = (dy) * inv_r3;
                force.x += Fx;
                force.y += Fy;
            }
        }
        return force;
    }
    int A = (coord.x >= node->x + node->size_x / 2); //find bounding box number
    int B = (coord.y >= node->y + node->size_y / 2);
    int C = 2*B + A;
    for(int i = 0; i < 4; ++i){
        if(i != C && node->children[i]){
            float x = (node->children[i])->COM.x;
            float y = (node->children[i])->COM.y;

            float dx = x - coord.x;
            float dy = y - coord.y;
            float r2 = dx*dx + dy*dy + softening;//1e-4f; // softening
            float r = sqrt(r2);
            float inv_r3 = 1.0f / (r2 * r);

            float Fx = ((node->children[i])->mass) * (x-coord.x) * inv_r3;
            float Fy = ((node->children[i])->mass)*(y-coord.y) * inv_r3;
            force.x += Fx;
            force.y += Fy;
        }
    }
    Vec2 subforce = GetForce(coord, node->children[C]);
    force.x += subforce.x;
    force.y += subforce.y;
    return force;
}

Vec2 Simulation::getForce(Vec2& coord, QuadNode* node){
    Vec2 force = {0.0f, 0.0f};
    if(node == nullptr){
        return {0.0, 0.0};
    }
    if(node->isleaf()){
        for(const Vec2& pos : node->same_node){
            float x = pos.x;
            float y = pos.y;

            float dx = x - coord.x;
            float dy = y - coord.y;
            if(dx*dx + dy*dy > 1e-6f){
                float r2 = dx*dx + dy*dy + softening;//1e-4f; // softening
                float r = sqrt(r2);
                float inv_r3 = 1.0f / (r2 * r);

                float Fx = (dx) * inv_r3;
                float Fy = (dy) * inv_r3;
                force.x += Fx;
                force.y += Fy;
            }
        }
        //cout << node->same_node.size() << '\n';
        //return force;
    }
    for(int i = 0; i < 4; ++i){
        if(node->children[i] != nullptr){
        float x = (node->children[i])->COM.x;
        float y = (node->children[i])->COM.y;
        
        float dx = x - coord.x;
        float dy = y - coord.y;
        float r2 = dx*dx + dy*dy + softening;//1e-4f; // softening
        if(((node->children[i]->size_x)*(node->children[i]->size_x)/r2 <= theta_squared) || node->children[i]->isleaf()){
            
            float r = sqrt(r2);
            float inv_r3 = 1.0f / (r2 * r);

            float Fx = ((node->children[i])->mass) * (x-coord.x) * inv_r3;
            float Fy = ((node->children[i])->mass)*(y-coord.y) * inv_r3;
            force.x += Fx;
            force.y += Fy;
        }
        else{
            //cout << "in" << '\n';
            Vec2 subforce = getForce(coord, node->children[i]);
            force.x += subforce.x;
            force.y += subforce.y;
        }
        }
    }
    //cout << '\n';
    return force;
}

void Simulation::UpdatePosition(Vec2& coord, Vec2& Force, Vec2& Velocity, Vec2& target){
    target.x = coord.x + ((Velocity.x * timestep) + (0.5 * Force.x * timestep * timestep));
    target.y = coord.y + ((Velocity.y * timestep) + (0.5 * Force.y * timestep * timestep));
    return; 
}

void Simulation::UpdateVelocity(Vec2& Force, Vec2& Velocity, Vec2& target){
    target.x = Velocity.x + (Force.x * timestep * 0.5);
    target.y = Velocity.y + (Force.y * timestep * 0.5);
    return; 
}

vector<Vec2> Simulation::generate_positions(){
    //std::mt19937 rng(42);  // fixed seed = reproducible
    std::mt19937 rng(std::random_device{}()); // random on every execution
    //std::uniform_real_distribution<float> dx(0.0f, windowDimensions.x);
    //std::uniform_real_distribution<float> dy(0.0f, windowDimensions.y);
    std::uniform_real_distribution<float> dx(700, 1300);
    std::uniform_real_distribution<float> dy(350, 650);
    std::vector<Vec2> pos;
    pos.reserve(N);

    for (int i = 0; i < N; ++i) {
        pos.push_back({ dx(rng), dy(rng) });
    }

    return pos;
}

vector<vector<Vec2>> Simulation::generateGalaxy(Vec2 core1, Vec2 core2, float M1, float M2){
    int R = 200;
    //std::mt19937 rng(42);  // fixed seed = reproducible
    std::mt19937 rng(std::random_device{}()); // random on every execution
    std::uniform_real_distribution<float> dr(10, R);
    std::uniform_real_distribution<float> da(0, 6.283);
    std::vector<Vec2> pos;
    std::vector<Vec2> vel;
    pos.reserve(N);
    vel.reserve(N);
    pos.push_back(core1);
    pos.push_back(core2);
    vel.push_back({0.0, 0.0});
    vel.push_back({0.0, 0.0});
    int val = 1;
    for (int i = 0; i < (N - 2) / 2; ++i) {
        float r = dr(rng);
        float a = da(rng);

        pos.push_back({
            r * cos(a) + core1.x,
            r * sin(a) + core1.y
        });
        vel.push_back({
            -1*sqrt(M1/sqrt(r*r + softening)) * sin(a),
            sqrt(M1/sqrt(r*r + softening)) * cos(a)
        });
    }

    for (int i = (N - 2) / 2; i < N - 2; ++i) {
        float r = dr(rng);
        float a = da(rng);

        pos.push_back({
            r * cos(a) + core2.x,
            r * sin(a) + core2.y
        });
        vel.push_back({
           sqrt(M2/sqrt(r*r + softening)) * sin(a),
            -1*sqrt(M2/sqrt(r*r + softening)) * cos(a)
        });
    }
    vector<vector<Vec2>> Config = {pos, vel};
    return Config;
}

void Simulation::initialise(vector<vector<Vec2>>& Config, vector<float>& mass){
    auto start = std::chrono::steady_clock::now();
    QuadNodePool nodes(2*N);
    QuadNode* root = nodes.allocate();
    root->x = corner.x;
    root->y = corner.y;
    root->size_x = windowDimensions.x;
    root->size_y = windowDimensions.y;
    for(int i = 0; i < N; ++i){
        build_tree(Config[0][i], root, nodes, mass[i]);
    }
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    std::cout << "Execution time for tree build: " << elapsed.count() << " ms\n";
    start = std::chrono::steady_clock::now();
    omp_set_num_threads(16);
    #pragma omp parallel for
    for(int i = 0; i < N; ++i){
        acc[i] = getForce(Config[0][i], root);
    }
    //nodes.freeList = nodes.block;
    end = std::chrono::steady_clock::now();
    elapsed = end - start;
    std::cout << "Execution time for force calc: " << elapsed.count() << " ms\n";
}

vector<vector<Vec2>> Simulation::NewState(vector<vector<Vec2>>& Config, vector<float>& mass, unsigned int& out){
    vector<vector<Vec2>> update;
    unsigned int nodecount = 0;
    auto start = std::chrono::steady_clock::now();
    QuadNodePool nodes(2*N);
    //nodes.freeList = nodes.block;
    QuadNode* root = nodes.allocate();
    root->x = corner.x;
    root->y = corner.y;
    root->size_x = windowDimensions.x;
    root->size_y = windowDimensions.y;
    for(int i = 0; i < N; ++i){
        UpdatePosition(Config[0][i], acc[i], Config[1][i], new_pos[i]);
        nodecount += build_tree(new_pos[i], root, nodes, mass[i]);
    }
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    std::cout << "Execution time for tree build: " << elapsed.count() << " ms\n";
    start = std::chrono::steady_clock::now();
    omp_set_num_threads(20);
    #pragma omp parallel for
    for(int i = 0; i < N; ++i){
        Vec2 temp;
        Vec2 force;
        force = getForce(new_pos[i], root);
        temp.x = force.x + acc[i].x;
        temp.y = force.y + acc[i].y;
        UpdateVelocity(temp, Config[1][i], new_vel[i]);
        acc[i].x = force.x;
        acc[i].y = force.y;
    }
    end = std::chrono::steady_clock::now();
    elapsed = end - start;
    std::cout << "Execution time for force: " << elapsed.count() << " ms\n";
    update.push_back(new_pos);
    update.push_back(new_vel);
    out = nodecount;
    return update;
}


