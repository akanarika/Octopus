// Octopus in C++
#include <GL/glew.h>
#include <GLUT/glut.h>
#include <algorithm>
#include <vector>
#include <glm/glm.hpp>
#include <sys/time.h>

// Phyxel Class
#include "phyxel.h"
// Class to load the object
// Reference from http://www.opengl-tutorial.org/
// beginners-tutorials/tutorial-7-model-loading/
#include "objloader.hpp"

/* World : seperate the space into small squares
* y
* ^    z
* |   /
* |  /
* | /
* | - - - - -> x
* 0
*/ 

// Window size
#define width 1024
#define height 800
#define grid_size 12
#define point_size 5

using namespace std;

//  time_step
float time_step =  0.25f;
float curr_time = 0;
double total_time = time_step;
struct timeval t1, t2;
double frame_timestep = 0;
float frame_time = 0;
float start_time = 0, fps = 0;
int frame_count = 0;

//  Selected vertex index
int selected_index = -1;

//  For camera control
float rot_x = 15;
float rot_y = 0;
float t_y = -3;
float dist = -23;
int mmb_click = 1;
glm::vec3 Up = glm::vec3(0,1,0);
glm::vec3 Right;
glm::vec3 view_dir;
GLint viewport[4];
GLdouble MV[16];
GLdouble P[16];

//  Define gravity vector
glm::vec3 gravity=glm::vec3(0.0f,-9.81f,0.0f);

//  Initialize for math
//  Poisson ratio
//  The Poisson's ratio of a stable, isotropic, 
//  linear elastic material will be 
//  greater than −1.0 or less than 0.5 because of 
//  the requirement for Young's modulus
float poisson =  0.4f;
//  Young's modulus
float young = 300.0f;
//  Material density
float density = 10000.f;
//  Spring constant   ( Energy = 1/2 kx^2, x = delta_distance)
float kv = 100;
//  Damping of velocity
float damping = 50.0f;

// Calculation for C
float d15 = young / (1.0f + poisson) / (1.0f - 2 * poisson);
float d16 = (1.0f - poisson) * d15;
float d17 = poisson * d15;
float d18 = young / 2 / (1.0f + poisson);
glm::vec3 D(d16, d17, d18);
glm::mat3 I = glm::mat3(1);

// Phyxels
vector<Phyxel*> phyxels;

// Loading model info
vector< glm::vec3 > vertices;
vector< glm::vec2 > uvs;
vector< glm::vec3 > normals;
bool res = loadOBJ("oct.obj", vertices, uvs, normals);

const int num_particle = vertices.size();
float scale = 0;

// Functions for GL
void OnRender();
void OnReshape(int nw, int nh);
void OnIdle();
void OnMouseMove(int x, int y);
void OnMouseDown(int button, int s, int x, int y);
void DrawGrid();
void InitGL();

// Functions for computing initial values
void GetNeighbors(int index, int k, vector<float>& dis, vector<Neighbor>& n);
void ComputeRadius(vector<float> &dists, float& r, float& h);
void SetNeighbors(int index);
void GetFactor();
void ComputeM(float dm, int index);
void ComputeRhoVolume(int index);
void ComputeInverseA(int index);

// Functions for compute update values
void ComputeLoop(float dt);
void UpdateF();
void ComputeJ();
void ComputeF();
void UpdatePos(float dt);

int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
    glutInitWindowSize(width, height);
    glutCreateWindow("Octopus");
    
    glutDisplayFunc(OnRender);
    glutReshapeFunc(OnReshape);
    glutIdleFunc(OnIdle);
    
    glutMouseFunc(OnMouseDown);
    glutMotionFunc(OnMouseMove);
    
    glewInit();
    
    // Initialize
    InitGL();
    
    glutMainLoop();
    
    return 0;
}

void OnRender() {
    int i = 0;
    float new_time = (float) glutGet(GLUT_ELAPSED_TIME);
    frame_time = new_time - curr_time;
    curr_time = new_time;
    gettimeofday(&t2, NULL);
    frame_timestep = (t2.tv_sec - t1.tv_sec) * 1000000.0;
    t1 = t2;
    total_time += frame_timestep;

    ++frame_count;
    
    if((new_time-start_time) > 100)
    {
        float elapsedTime = new_time - start_time;
        fps = (frame_count / elapsedTime) * 100000.0;
        start_time = new_time;
        frame_count = 0;
    }
    
    glutSetWindowTitle("Octopus");
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    glTranslatef(0, t_y, dist);
    glRotatef(rot_x, 1, 0, 0);
    glRotatef(rot_y, 0, 1, 0);
    
    glGetDoublev(GL_MODELVIEW_MATRIX, MV);
    view_dir.x = (float)-MV[2];
    view_dir.y = (float)-MV[6];
    view_dir.z = (float)-MV[10];
    Right = glm::cross(view_dir, Up);
    
    //draw grid
    DrawGrid();

    //draw points
    glBegin(GL_POINTS);
    for(i = 0; i < num_particle; i++) {
        glm::vec3 p = phyxels[i]->X;
        
        if(i == selected_index)
            glColor3f(1,1,1);
        else
            glColor3f(0.3,0.3,0.3);
            
        glVertex3f(p.x,p.y,p.z);
    }
    
    glEnd();
    
    glutSwapBuffers();
}

void DrawGrid()
{
    glBegin(GL_LINES);
    glColor3f(0.05f, 0.05f, 0.05f);
    for(int i = -grid_size; i <= grid_size; i++)
    {
        glVertex3f((float)i, 0, (float)-grid_size);
        glVertex3f((float)i, 0, (float)grid_size);
        
        glVertex3f((float)-grid_size, 0, (float)i);
        glVertex3f((float)grid_size, 0, (float)i);
    }
    glEnd();
}

void OnReshape(int nw, int nh) {
    glViewport(0,0,nw, nh);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60, (GLfloat)nw / (GLfloat)nh, 1.f, 100.0f);
    
    glGetIntegerv(GL_VIEWPORT, viewport);
    glGetDoublev(GL_PROJECTION_MATRIX, P);
    
    glMatrixMode(GL_MODELVIEW);
}

void OnIdle() {
    
    if ( total_time >= time_step )
    {
        ComputeLoop(time_step);
        total_time -= time_step;
    }
    
    glutPostRedisplay();
}

int prev_x=0, prev_y=0;
void OnMouseDown(int button, int s, int x, int y)
{
    if (s == GLUT_DOWN)
    {
        prev_x = x;
        prev_y = y;

        GLfloat window_x, window_y, window_z;
        GLdouble obj_x, obj_y, obj_z;

        window_x = (float) x;
        window_y = viewport[3] - (float) y;

        glReadPixels(x, (int)window_y, 1, 1, GL_DEPTH_COMPONENT, 
                     GL_FLOAT, &window_z);
        
        gluUnProject(window_x, window_y, window_z, MV, P, 
                     viewport, &obj_x, &obj_y, &obj_z);
        glm::vec3 pt(obj_x, obj_y, obj_z);
        int i = 0;
        for(i = 0; i < num_particle; i++) {
            if(glm::distance(phyxels[i]->X,pt) < 0.1) {
                selected_index = i;
                printf("Intersected at %d\n",i);
                break;
            }
        }
    }
    
    if(button == GLUT_MIDDLE_BUTTON)
        mmb_click = 0;
    else
        mmb_click = 1;
    
    if(s==GLUT_UP) {
        selected_index= -1;
        glutSetCursor(GLUT_CURSOR_INHERIT);
    }
}

void OnMouseMove(int x, int y)
{
    if(selected_index == -1) {
        if (mmb_click == 0)
            dist *= (1 + (y - prev_y)/60.0f);
        else
        {
            rot_y += (x - prev_x)/5.0f;
            rot_x += (y - prev_y)/5.0f;
        }
    } else {
        float delta = 1500/abs(dist);
        float valX = (x - prev_x)/delta;
        float valY = (prev_y - y)/delta;
        if(abs(valX)>abs(valY))
            glutSetCursor(GLUT_CURSOR_LEFT_RIGHT);
        else
            glutSetCursor(GLUT_CURSOR_UP_DOWN);
        
        phyxels[selected_index]->v = glm::vec3(0);
        phyxels[selected_index]->X += Right[0]*valX ;
        float newValue = phyxels[selected_index]->X.y+Up[1]*valY;
        if(newValue>0)
            phyxels[selected_index]->X.y = newValue;
        phyxels[selected_index]->X.z += Right[2]*valX + Up[2]*valY;
        
        phyxels[selected_index]->v = glm::vec3(0);
    }
    prev_x = x;
    prev_y = y;
    
    glutPostRedisplay();
}

void InitGL() {
    start_time = (float)glutGet(GLUT_ELAPSED_TIME);
    curr_time = start_time;
    
    gettimeofday(&t1, NULL);
    
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_POINT_SMOOTH);
    int i = 0, j = 0, k = 0, count=0;
    
    float ypos = 4.0f;

    for (i = 0; i < num_particle; i++) {
        Phyxel* phyxel = new Phyxel();
        phyxels.push_back(phyxel);
    }
    
    printf("vertices size: %lu\n", vertices.size());

    for (i = 0; i < vertices.size(); i++) {
        phyxels[i]->setX(glm::vec3(vertices[i].x * 3,
                    vertices[i].y * 3 + 4, vertices[i].z * 3));
        phyxels[i]->setXi(phyxels[i]->X);
        // phyxels[i]->isFixed = (vertices[i].y > 2);
    }
    
    
    for (i = 0; i < num_particle; ++i) {
        vector<float> dist;
        GetNeighbors(i, 10, dist, phyxels[i]->neighbors);
        ComputeRadius(dist, phyxels[i]->r, phyxels[i]->h);
    }
    
    for (i = 0; i < num_particle; i++) {
        SetNeighbors(i);
    }
    
    GetFactor();
    
    for(i=0; i < num_particle; ++i)
        ComputeM(density, i);
    
    for(i=0; i < num_particle; ++i)
        ComputeRhoVolume(i);
    
    for(i=0; i < num_particle; ++i) {
        ComputeInverseA(i);
    }

    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glPointSize(point_size );
    
}

typedef std::pair<int, float> idx_to_dist;

struct Cmp {
    bool operator()(const idx_to_dist &lhs, const idx_to_dist &rhs) {
        return lhs.second < rhs.second;
    }
};

void GetNeighbors(int index, int k, vector<float>& dis, vector<Neighbor>& n) {
    vector<idx_to_dist> distances;

    for(int i = 0; i < num_particle; i++) {
        if(index!=i) {
            idx_to_dist m;
            m.first = i;
            m.second= fabs(glm::distance(phyxels[index]->X, phyxels[i]->X));
            distances.push_back(m);
        }
    }

    sort (distances.begin(), distances.end(), Cmp());

    for(int i = 0; i < k; i++) {
        Neighbor new_neighbor;
        new_neighbor.j = distances[i].first;
        dis.push_back(distances[i].second);
        n.push_back(new_neighbor);
    }
}

void ComputeRadius(vector<float> &dists, float& r, float& h)
{
    float avg = 0.f;
    for(size_t i=0;i<dists.size();i++)
        avg += dists[i];
    r = avg / dists.size();

    h = 3.0f * r;
}

void SetNeighbors(int index)
{
    float mul_factor = float(315.0f
                            /(64*M_PI*pow(phyxels[index]->h,9.0f)));
    float h2 = phyxels[index]->h * phyxels[index]->h;
    for(size_t i=0;i<phyxels[index]->neighbors.size();i++)
    {
        Neighbor& n = phyxels[index]->neighbors[i];
        n.rdist  = phyxels[n.j]->Xi - phyxels[index]->Xi;
        float r2 = glm::dot(n.rdist, n.rdist);

        n.w = mul_factor * pow(h2 - r2, 3);
    }
}

void GetFactor() {
    scale = 0.f;
    int i=0;
    for ( i =0; i < num_particle; i++ ) {
        float sum = 0.f;

        vector<Neighbor>& pNeighbor = phyxels[i]->neighbors;
        for(size_t j = 0; j < pNeighbor.size(); j++) {
            sum += pow(phyxels[pNeighbor[j].j]->r, 3) * pNeighbor[j].w;
        }
        scale += 1.0f / sum;
    }
    scale /= num_particle;
}

void ComputeM(float dm, int index)
{
    phyxels[index]->m = scale * pow(phyxels[index]->r, 3) * dm;
}

void ComputeRhoVolume(int index)
{
    phyxels[index]->rho = 0.f;
    vector<Neighbor>& pNeighbor = phyxels[index]->neighbors;
    for(size_t i = 0; i < pNeighbor.size(); i++)
        phyxels[index]->rho += phyxels[pNeighbor[i].j]->m * pNeighbor[i].w;       // Eq. 10 page 4
    phyxels[index]->vol =  phyxels[index]->m / phyxels[index]->rho;
}

// A_i = sum_j(x_ij * x_ij^T * w_ij)
void ComputeInverseA(int index)
{
    glm::mat3 A, A_sum, V;
    A = glm::mat3(0);

    vector<Neighbor>& pNeighbor = phyxels[index]->neighbors;
    for(size_t i = 0;i < pNeighbor.size();i++)
    {
        A_sum = glm::outerProduct(pNeighbor[i].rdist, 
                pNeighbor[i].rdist * pNeighbor[i].w);
        A += A_sum;
    }

    if(glm::determinant(A) != 0.0)
        phyxels[index]->Ai = glm::inverse(A);

    phyxels[index]->di = glm::vec3(0);

    for(size_t i = 0; i < pNeighbor.size(); i++)
    {
        glm::vec3 Xij_Wij = pNeighbor[i].rdist * pNeighbor[i].w;
        pNeighbor[i].dj = phyxels[index]->Ai * Xij_Wij;
        phyxels[index]->di -= pNeighbor[i].dj;
    }
}

void ComputeLoop(float dt) {
    UpdateF();
    UpdatePos(dt);
}

void UpdateF()
{
    for(int i = 0; i < num_particle; i++) {
        if(!phyxels[i]->isFixed)
            phyxels[i]->F = gravity;
        else
            phyxels[i]->F = glm::vec3(0);

        phyxels[i]->F -= phyxels[i]->v * damping;
    }

    ComputeJ();

    ComputeF();

    for(int i = 0; i < num_particle; i++) {
        glm::mat3 F_e, F_v;
        F_e =  -2 * phyxels[i]->vol * (phyxels[i]->J * phyxels[i]->sigma) ;
        glm::vec3 J_u = glm::vec3(phyxels[i]->J[0][0], 
                                  phyxels[i]->J[0][1],
                                  phyxels[i]->J[0][2]);
        glm::vec3 J_v = glm::vec3(phyxels[i]->J[1][0],
                                  phyxels[i]->J[1][1],
                                  phyxels[i]->J[1][2]);
        glm::vec3 J_w = glm::vec3(phyxels[i]->J[2][0],
                                  phyxels[i]->J[2][1],
                                  phyxels[i]->J[2][2]);

        glm::vec3 row0 = glm::cross(J_v, J_w);
        glm::vec3 row1 = glm::cross(J_w, J_u);
        glm::vec3 row2 = glm::cross(J_u, J_v);

        glm::mat3 M = glm::transpose(glm::mat3(row0, row1, row2));

        F_v = -phyxels[i]->vol * kv 
              * (glm::determinant(phyxels[i]->J) - 1) * M ;

        vector<Neighbor>& pNeighbor = phyxels[i]->neighbors;
        for(size_t j = 0; j < pNeighbor.size(); j++)
            phyxels[pNeighbor[j].j]->F += (F_e + F_v) * pNeighbor[j].dj;

        phyxels[i]->F += (F_e + F_v) * phyxels[i]->di;
    }
}

void ComputeJ()
{
    for(int i = 0; i < num_particle; i++) {
        vector<Neighbor>& pNeighbor = phyxels[i]->neighbors;
        for(size_t j = 0; j < pNeighbor.size(); j++)
            phyxels[pNeighbor[j].j]->U = phyxels[pNeighbor[j].j]->X 
                                        - phyxels[pNeighbor[j].j]->Xi;

        glm::mat3 B = glm::mat3(0);

        glm::mat3 du = glm::mat3(0);
        glm::mat3 du_tr = glm::mat3(0);

        for(size_t j = 0; j < pNeighbor.size(); j++)
        {
            glm::mat3 bj = glm::mat3(0);
            bj = glm::outerProduct(phyxels[pNeighbor[j].j]->U -
                                       phyxels[i]->U, 
                                       pNeighbor[j].rdist * pNeighbor[j].w );
            B += bj;
        }
        B = glm::transpose(B);

        du = phyxels[i]->Ai * B;
        du_tr = glm::transpose(du);
        phyxels[i]->J=glm::mat3(1);
        phyxels[i]->J += du_tr;
    }

}

void ComputeF()
{
    for(int i = 0; i < num_particle; i++) {
        glm::mat3 Jtr = glm::transpose(phyxels[i]->J);
        phyxels[i]->epsilon = (Jtr * phyxels[i]->J) - I;

        glm::mat3& e = phyxels[i]->epsilon;  // epsilon
        glm::mat3& s = phyxels[i]->sigma;  // sigma

        s[0][0] = D.x * e[0][0] 
                  + D.y * e[1][1] 
                  + D.y * e[2][2];
        s[1][1] = D.y * e[0][0]
                  + D.x * e[1][1]
                  + D.y*e[2][2];
        s[2][2] = D.y*e[0][0]
                  + D.y*e[1][1] 
                  + D.x*e[2][2];

        s[0][1] = D.z * e[0][1];
        s[1][2] = D.z * e[1][2];
        s[2][0] = D.z * e[2][0];

        s[0][2] = s[2][0];
        s[1][0] = s[0][1];
        s[2][1] = s[1][2];
    }
}

void UpdatePos(float dt) {
    float dt2 = dt * dt;
    float half_dt2 = dt2 * 0.5f;
    
    int i = 0;
    for(i = 0; i < num_particle; i++) {
        if(!phyxels[i]->isFixed) {
            phyxels[i]->acc = phyxels[i]->F / phyxels[i]->m;

            phyxels[i]->X += dt * phyxels[i]->v 
                             + (phyxels[i]->acc * half_dt2);
            
            if(phyxels[i]->X.y < 0) {
                phyxels[i]->X.y = 0;
            }
        }
    }
    
    UpdateF();

    for(i=0;i<num_particle;i++) {
        if(!phyxels[i]->isFixed)
            phyxels[i]->v += ((phyxels[i]->F/phyxels[i]->m 
                + phyxels[i]->acc)*half_dt2);
    }

}