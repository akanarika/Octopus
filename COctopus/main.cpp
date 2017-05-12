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

// Window size
#define width 1024
#define height 800
#define grid_size 12
#define point_size 5

using namespace std;

// time_step
float time_step =  0.25f;
float curr_time = 0;
double total_time = time_step;

// Selected vertex index
int selected_index = -1;

// For camera control
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

// Define gravity vector
glm::vec3 gravity=glm::vec3(0.0f,-9.81f,0.0f);

// Time step
struct timeval t1, t2;
double frame_timestep = 0;
float frame_time = 0;
float start_time = 0, fps = 0;
int frame_count = 0;

// Initialize for math
// Poisson ratio
float poisson =  0.4f;
// Young modulus
float young = 300.0f;
// Material density
float density = 10000.f;
// K
float kv = 100;
// Damping of velocity
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
        
        if(i==selected_index)
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
        glVertex3f((float)i,0,(float)-grid_size);
        glVertex3f((float)i,0,(float)grid_size);
        
        glVertex3f((float)-grid_size,0,(float)i);
        glVertex3f((float)grid_size,0,(float)i);
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
        GLdouble objX, objY, objZ;

        window_x = (float) x;
        window_y = viewport[3] - (float) y;

        glReadPixels(x, (int)window_y, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &window_z);
        
        gluUnProject(window_x, window_y, window_z, MV, P, viewport, &objX, &objY, &objZ);
        glm::vec3 pt(objX, objY, objZ);
        int i=0;
        for(i=0;i<num_particle;i++) {
            if( glm::distance(phyxels[i]->X,pt)<.1) {
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

    //Get the distances of current point to all other points
    for(int i = 0; i < num_particle; i++) {
        if(index!=i) {
            idx_to_dist m;
            m.first = i;
            m.second= fabs(glm::distance(phyxels[index]->X, phyxels[i]->X));
            distances.push_back(m);
        }
    }

    //sort the distances
    sort (distances.begin(), distances.end(), Cmp());

    //now take the top k neighbors based on distance
    for(int i = 0; i < k; i++) {
        Neighbor ne;
        ne.j = distances[i].first;
        dis.push_back(distances[i].second);
        n.push_back(ne);
    }
}

void ComputeRadius(vector<float> &dists, float& r, float& h)
{
    //For this read section 3.2 Initialization on page 4
    //look through all neigbour distances
    //and average the distance. This will give us r
    float avg = 0.f;
    for(size_t i=0;i<dists.size();i++)
        avg += dists[i];
    r = avg / dists.size();

    // compute the support radii H = 3 * R
    h = 3.0f * r;
}

void SetNeighbors(int index)
{
    //For this read section 3.2 Initialization on page 4
    //Based on Eq. 9
    float mul_factor = float(315.0f
                            /(64*M_PI*pow(phyxels[index]->h,9.0f)));//35 / (float)(32 * M_PI * pow(h, 7));
    float h2 = phyxels[index]->h * phyxels[index]->h;
    for(size_t i=0;i<phyxels[index]->neighbors.size();i++)
    {
        Neighbor& n = phyxels[index]->neighbors[i];
        n.rdist  = phyxels[n.j]->Xi - phyxels[index]->Xi;
        float r2 = glm::dot(n.rdist, n.rdist);      //r = sqrt(Xij.x*Xij.x+Xij.y*Xij.y+Xij.z*Xij.z)
                                                    //r^2 = dot(Xij,Xij);
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
    // This is the common scaling factor to compute the mass of each phyxel.
    // See last paragraph of Section 3.2 on page 4
    scale /= num_particle;
}

void ComputeM(float dm, int index)
{
    // See last paragraph of Section 3.2 on page 4
    phyxels[index]->m = scale * pow(phyxels[index]->r, 3) * dm;
}

void ComputeRhoVolume(int index)
{
    // See last paragraph of Section 3.2 on page 4
    phyxels[index]->rho = 0.f;
    vector<Neighbor>& pNeighbor = phyxels[index]->neighbors;
    for(size_t i = 0; i < pNeighbor.size(); i++)
        phyxels[index]->rho += phyxels[pNeighbor[i].j]->m * pNeighbor[i].w;       // Eq. 10 page 4
    phyxels[index]->vol =  phyxels[index]->m / phyxels[index]->rho;
}

void ComputeInverseA(int index)
{
    glm::mat3 A, A_sum, V;
    A =glm::mat3(0);

    vector<Neighbor>& pNeighbor = phyxels[index]->neighbors;
    for(size_t i=0;i < pNeighbor.size();i++)
    {
        A_sum = glm::outerProduct(pNeighbor[i].rdist, 
                pNeighbor[i].rdist * pNeighbor[i].w); // Eq. 14
        A += A_sum;
    }

    if(glm::determinant(A) != 0.0)
        phyxels[index]->Ai = glm::inverse(A);      // Eq. 14, inverted moment matrix
    else
    {
        // if MomentMatrix is not invertible it means that there are less than 4 neighbors
        // or the neighbor phyxels are coplanar or colinear
        // We should use SVD to extract the inverse but I have left this as is.
        //Read section 3.3 last paragraph
        puts("Warning: Singular matrix!!!");
    }

    phyxels[index]->di=glm::vec3(0);

    for(size_t i=0;i<pNeighbor.size();i++)
    {
        glm::vec3 Xij_Wij = pNeighbor[i].rdist * pNeighbor[i].w;
        pNeighbor[i].dj = phyxels[index]->Ai * Xij_Wij;    //Eq. 21 page 5
        phyxels[index]->di -= pNeighbor[i].dj;               //Eq. 20 page 5
    }
}

void ComputeLoop(float dt) {
    UpdateF();
    UpdatePos(dt);
}

void UpdateF()
{
    //Calculate external force
    for(int i=0;i<num_particle;i++) {
        if(!phyxels[i]->isFixed)
            phyxels[i]->F = gravity;
        else
            phyxels[i]->F = glm::vec3(0);

        //Add velocity damping
        phyxels[i]->F -= phyxels[i]->v * damping;
    }

    ComputeJ();

    ComputeF();


    //Calculate internal force using the stress and Jacobians
    for(int i=0;i<num_particle;i++) {
        glm::mat3 F_e, F_v;
        F_e =  -2 * phyxels[i]->vol * (phyxels[i]->J * phyxels[i]->sigma) ;    // Eq. 18
        glm::vec3 J_u = glm::vec3(phyxels[i]->J[0][0], phyxels[i]->J[0][1],phyxels[i]->J[0][2]);
        glm::vec3 J_v = glm::vec3(phyxels[i]->J[1][0], phyxels[i]->J[1][1],phyxels[i]->J[1][2]);
        glm::vec3 J_w = glm::vec3(phyxels[i]->J[2][0], phyxels[i]->J[2][1],phyxels[i]->J[2][2]);

        glm::vec3 row0 = glm::cross(J_v, J_w);  //Eq. 22
        glm::vec3 row1 = glm::cross(J_w, J_u);  //Eq. 22
        glm::vec3 row2 = glm::cross(J_u, J_v);  //Eq. 22
        glm::mat3 M = glm::transpose(glm::mat3(row0, row1, row2));   //Eq. 22

        F_v = -phyxels[i]->vol * kv * (glm::determinant(phyxels[i]->J) - 1) * M ; //Eq. 22

        vector<Neighbor>& pNeighbor = phyxels[i]->neighbors;
        for(size_t j=0;j<pNeighbor.size();j++)
            phyxels[pNeighbor[j].j]->F += (F_e + F_v) * pNeighbor[j].dj;

        phyxels[i]->F += (F_e + F_v) * phyxels[i]->di;
    }
}

void ComputeJ()
{
    for(int i=0;i<num_particle;i++) {
        vector<Neighbor>& pNeighbor = phyxels[i]->neighbors;
        for(size_t j=0;j<pNeighbor.size();j++)
            phyxels[pNeighbor[j].j]->U = phyxels[pNeighbor[j].j]->X 
                                        - phyxels[pNeighbor[j].j]->Xi;

        glm::mat3 B=glm::mat3(0);       // help matrix used to compute the sum in Eq. 15

        // reset du and du_tr
        glm::mat3 du=glm::mat3(0);
        glm::mat3 du_tr=glm::mat3(0);

        for(size_t j=0;j < pNeighbor.size();j++)
        {
            glm::mat3 Bj=glm::mat3(0);
            //Eq. 15 right hand side terms with A_inv
            Bj =glm::outerProduct(phyxels[pNeighbor[j].j]->U - phyxels[i]->U, 
                            pNeighbor[j].rdist * pNeighbor[j].w );
            B += Bj;
        }
        B = glm::transpose(B);

        du = phyxels[i]->Ai * B;   // Eq. 15 page 4
        du_tr = glm::transpose(du);
        phyxels[i]->J=glm::mat3(1);
        phyxels[i]->J += du_tr;      // Eq. 1
    }

}

void ComputeF()
{
    for(int i=0;i<num_particle;i++) {
        glm::mat3 Jtr = glm::transpose(phyxels[i]->J);
        phyxels[i]->epsilon = (Jtr * phyxels[i]->J) - I;      // formula 3, Green-Saint Venant non-linear tensor

        glm::mat3& e= phyxels[i]->epsilon;
        glm::mat3& s= phyxels[i]->sigma;

        s[0][0] = D.x*e[0][0]+D.y*e[1][1]+D.y*e[2][2];
        s[1][1] = D.y*e[0][0]+D.x*e[1][1]+D.y*e[2][2];
        s[2][2] = D.y*e[0][0]+D.y*e[1][1]+D.x*e[2][2];

        s[0][1] = D.z*e[0][1];
        s[1][2] = D.z*e[1][2];
        s[2][0] = D.z*e[2][0];

        s[0][2] = s[2][0];
        s[1][0] = s[0][1];
        s[2][1] = s[1][2];
    }
}

void UpdatePos(float dt) {
    int i=0;
    float dt2 = dt*dt;
    float half_dt2 = dt2*0.5f;

    for(i=0;i<num_particle;i++) {

        //X_(i+1) = X_i + V_i*dt + 1/2*a_i*dt^2
        if(!phyxels[i]->isFixed) {
            phyxels[i]->acc = phyxels[i]->F/phyxels[i]->m;

            phyxels[i]->X += dt*phyxels[i]->v+(phyxels[i]->acc*half_dt2);
            //printf("acc %d: pos (%f,%f,%f)\n", i, acc0[i].x, acc0[i].y, acc0[i].z);

            if(phyxels[i]->X.y <0) {
                phyxels[i]->X.y = 0;
            }
            //printf("X %d: pos (%f,%f,%f)", i, X[i].x, X[i].y, X[i].z);
        }
    }
    //Calculate the new acceleration
    UpdateF();

    //V_(i+1) = V_i + ((a_i+a_(i+1))/2)*dt^2
    for(i=0;i<num_particle;i++) {
        if(!phyxels[i]->isFixed)
            phyxels[i]->v += ((phyxels[i]->F/phyxels[i]->m 
                + phyxels[i]->acc)*half_dt2);
    }

}