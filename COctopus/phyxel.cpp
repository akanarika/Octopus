#include "phyxel.h"

Phyxel::Phyxel() {
    Xi = glm::vec3(0);
    X = glm::vec3(0);
    U = glm::vec3(0);
    v = glm::vec3(0);
    isFixed = false;

    Ai = glm::mat3(0);

    di = glm::vec3(0);
    sigma = glm::mat3(0);
    epsilon = glm::mat3(0);
    J = glm::mat3(0);

    F = glm::vec3(0);
}

void Phyxel::setX(glm::vec3 x) {
    X = x;
}


void Phyxel::setXi(glm::vec3 xi) {
    Xi = xi;
}