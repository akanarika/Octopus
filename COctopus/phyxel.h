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
    glm::vec3 Xi; //init pos
    glm::vec3 X;  //curr pos
    glm::vec3 U;  //displacements

    float m; //mass
    glm::mat3 Ai;

    glm::vec3 di;
    glm::mat3 sigma;  //stress
    glm::mat3 epsilon;  //strains
    glm::mat3 J;  //Jacobian
    glm::vec3  acc; 

    vector<Neighbor> neighbors;

    float r;  //distance to neighbors
    float h;  //support radius
    float rho;  //density
    float vol;  //volume
    glm::vec3 v;  //current velocity
    glm::vec3 F;  //force
    bool isFixed;
};



#endif