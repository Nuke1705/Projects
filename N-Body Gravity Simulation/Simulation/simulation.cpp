//#include <SFML/Graphics.hpp>
#include <vector>
#include <random>
#include "simulation.h"

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
    vector<Vec2> same_node;
    QuadNode* nxtfree = nullptr;
    //node functions
    bool isleaf(){
        return (children[0]==nullptr)&&(children[1]==nullptr)&&(children[2]==nullptr)&&(children[3]==nullptr);
    }
    void set(float X, float Y, float size_X, float size_Y){
        x = X;
        y = Y;
        size_x = size_X;
        size_y = size_Y;
        mass = 0;
        nxtfree = nullptr;
        fill(children, children + 4, nullptr);
        return;
    }

};

class QuadNodePool{
public:
    QuadNodePool(size_t blockSize = 1024)
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

Simulation::Simulation(const Vec2& Dimensions, long long ParticleCount, float timestep)
: windowDimensions(Dimensions), N(ParticleCount), timestep(timestep){
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

unsigned int Simulation::build_tree(Vec2& coord, QuadNode* node, QuadNodePool& nodes){
    unsigned int nodecount = 0;
    if(node->mass == 0){
        node->mass = 1;   //root case
        node->COM = coord;
        nodecount++;
        return nodecount;
    }
    else{
        if(node->isleaf()){ //if the node has a particle but no children i.e m = 1
            if(min(node->size_x/2, node->size_y/2) < MIN_SIZE){
                node->same_node.push_back(coord);
                node->mass++;
                Vec2 ncom = com(coord, 1.0, node->COM, node->mass);
                node->COM = ncom; 
                return nodecount;
            }
            QuadNode* branch1 = nodes.allocate();
            branch1->COM = node->COM;
            int a = (node->COM.x >= node->x + node->size_x / 2);
            int b = (node->COM.y >= node->y + node->size_y / 2);
            int c = 2*b + a;
            node->children[c] = branch1;
            branch1->mass = 1;
            branch1->x = node->x + a*(node->size_x/2);
            branch1->y = node->y + b*(node->size_y/2);
            branch1->size_x = node->size_x/2;
            branch1->size_y = node->size_y/2;

            int A = (coord.x >= node->x + node->size_x / 2);
            int B = (coord.y >= node->y + node->size_y / 2);
            int C = 2*B + A;
            
            if(C == c){
               nodecount += build_tree(coord, branch1, nodes);
            }
            else{
                QuadNode* branch2 = nodes.allocate();
                branch2->COM = coord;
                node->children[C] = branch2;
                branch2->mass = 1;
                branch2->x = node->x + A*(node->size_x/2);
                branch2->y = node->y + B*(node->size_y/2);
                branch2->size_x = node->size_x/2;
                branch2->size_y = node->size_y/2;
                nodecount++;
            }
            
            
            Vec2 ncom = com(coord, 1.0, node->COM, node->mass);
            node->COM = ncom;
            node->mass++;
            nodecount++;
            return nodecount;
            
        }
        else{ //has atleast one child i.e. mass>1
            int A = (coord.x >= node->x + node->size_x / 2);
            int B = (coord.y >= node->y + node->size_y / 2);
            int C = 2*B + A;
            Vec2 ncom = com(coord, 1.0, node->COM, node->mass);
            node->COM = ncom;
            node->mass++;
            if(node->children[C] == nullptr){
                QuadNode* branch2 = nodes.allocate();
                branch2->COM = coord;
                node->children[C] = branch2;
                branch2->mass = 1;
                branch2->x = node->x + A*(node->size_x/2);
                branch2->y = node->y + B*(node->size_y/2);
                branch2->size_x = node->size_x/2;
                branch2->size_y = node->size_y/2;
                nodecount++;
            }
            else{
               nodecount += build_tree(coord, node->children[C], nodes);
            }
            return nodecount;
        }
    }

}

Vec2 Simulation::GetForce(Vec2& coord, QuadNode* node){
    int A = (coord.x >= node->x + node->size_x / 2); //find bounding box number
    int B = (coord.y >= node->y + node->size_y / 2);
    int C = 2*B + A;
    Vec2 force = {0.0f, 0.0f};
    for(int i = 0; i < 4; ++i){
        if(i != C && node->children[i]){
            float x = (node->children[i])->COM.x;
            float y = (node->children[i])->COM.y;

            float dx = x - coord.x;
            float dy = y - coord.y;
            float r2 = dx*dx + dy*dy + 25.0f;//1e-4f; // softening
            float r = sqrt(r2);
            float inv_r3 = 1.0f / (r2 * r);

            float Fx = ((node->children[i])->mass) * (x-coord.x) * inv_r3;
            float Fy = ((node->children[i])->mass)*(y-coord.y) * inv_r3;
            force.x += Fx;
            force.y += Fy;
        }
    }
    if(node->isleaf()){
        for(const Vec2& pos : node->same_node){
            float x = pos.x;
            float y = pos.y;

            float dx = x - coord.x;
            float dy = y - coord.y;
            if(dx*dx + dy*dy < 1e-6f){
                float r2 = dx*dx + dy*dy + 0.1f;//1e-4f; // softening
                float r = sqrt(r2);
                float inv_r3 = 1.0f / (r2 * r);

                float Fx = (x-coord.x) * inv_r3;
                float Fy = (y-coord.y) * inv_r3;
                force.x += Fx;
                force.y += Fy;
            }
        }
        return force;
    }
    else{
        Vec2 subforce = GetForce(coord, node->children[C]);
        force.x += subforce.x;
        force.y += subforce.y;
        return force;
    }
}

Vec2 Simulation::UpdatePosition(Vec2& coord, Vec2& Force, Vec2& Velocity){
    Vec2 new_coord = {0.0f, 0.0f};
    new_coord.x += coord.x + ((Velocity.x * timestep) + (0.5 * Force.x * timestep * timestep));
    new_coord.y += coord.y + ((Velocity.y * timestep) + (0.5 * Force.y * timestep * timestep));
    if(new_coord.x >= windowDimensions.x){
        new_coord.x = windowDimensions.x;
    }
    if(new_coord.y >= windowDimensions.y){
        new_coord.y = windowDimensions.y;
    }
    if(new_coord.x <= 0.0f){
        new_coord.x = 0.0f;
    }
    if(new_coord.y <= 0.0f){
        new_coord.y = 0.0f;
    }
    return new_coord; 
}

Vec2 Simulation::UpdateVelocity(Vec2& Force, Vec2& Velocity, Vec2& coord){
    Vec2 new_vel = {0.0f, 0.0f};
    if(coord.x == windowDimensions.x || coord.y == windowDimensions.y || coord.x == 0.0f || coord.y == 0.0f){
        new_vel.x = 0.0f;
        new_vel.y = 0.0f;
        return new_vel;
    }
    new_vel.x += Velocity.x + (Force.x * timestep);
    new_vel.y += Velocity.y + (Force.y * timestep);
    return new_vel; 
}



vector<Vec2> Simulation::generate_positions(){
    //std::mt19937 rng(42);  // fixed seed = reproducible
    std::mt19937 rng(std::random_device{}()); // random on every execution
    std::uniform_real_distribution<float> dx(0.0f, windowDimensions.x);
    std::uniform_real_distribution<float> dy(0.0f, windowDimensions.y);

    std::vector<Vec2> pos;
    pos.reserve(N);

    for (int i = 0; i < N; ++i) {
        pos.push_back({ dx(rng), dy(rng) });
    }

    return pos;
}

vector<vector<Vec2>> Simulation::NewState(vector<vector<Vec2>>& Config, unsigned int& out){
    vector<vector<Vec2>> update;
    unsigned int nodecount = 0;
    {
    vector<Vec2> new_pos(N);
    vector<Vec2> new_vel(N);
    QuadNodePool nodes;
    QuadNode* root = nodes.allocate();
    root->set(0.0f, 0.0f, windowDimensions.x, windowDimensions.y);
    for(Vec2& pos : Config[0]){
        nodecount += build_tree(pos, root, nodes);
    }
    for(int i = 0; i < N; ++i){
        Vec2 force = GetForce(Config[0][i], root);
        new_pos[i] = UpdatePosition(Config[0][i], force, Config[1][i]);
        new_vel[i] = UpdateVelocity(force, Config[1][i], Config[0][i]);
    }
    update.push_back(new_pos);
    update.push_back(new_vel);
    }
    out = nodecount;
    return update;
}


