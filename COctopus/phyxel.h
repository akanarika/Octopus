#include <vector>
#include <glm/glm.hpp>
#ifndef PHYXEL_H
#define PHYXEL_H

using namespace std;

struct Neighbor { 
    int  j;
    float w;    
    glm::vec3 rdist; 
    glm::vec3 dj;
};

class Phyxel {
public:
    Phyxel();
    void setX(glm::vec3 pos);
    void setXi(glm::vec3 pos);
    glm::vec3 Xi;  // Original material coordinates
    glm::vec3 X;  // Current coordinates
    glm::vec3 U;  // Displacement

    float m;  // Mass
    glm::mat3 Ai;  // Inversed A Matirx

    glm::vec3 di;
    glm::mat3 sigma;
    glm::mat3 epsilon;
    glm::mat3 J;  // Jacobian
    glm::vec3  acc;  // Acceleration

    vector<Neighbor> neighbors;

    float r;  // Distance to neighbors
    float h;  // Support radius
    float rho;  // Density
    float vol;  // Volume
    glm::vec3 v;  // Current velocity
    glm::vec3 F;  // Force
    bool isFixed;  // Whether it's fixed
};


#endif