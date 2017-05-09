// Octopus in C++
#include <GL/glew.h>
#include <GL/glut.h>
#include <algorithm>
#include <vector>
#include <glm/glm.hpp>
#include <sys/time.h>
#include "phyxel.h"

#define _USE_MATH_DEFINES

using namespace std;

const int width = 1024, height = 1024;

int numX = 5, numY=20, numZ=5;
const int num_particle = numX*numY*numZ;
int sizeX = 1;
int sizeY = 4;
int sizeZ = 1;
float hsizeX = sizeX/2.0f;
float hsizeY = sizeY/2.0f;
float hsizeZ = sizeZ/2.0f;

bool bShowForces=false, bShowJacobians=false;

float timeStep =  10/60.0f;
float currentTime = 0;
double accumulator = timeStep;
int selected_index = -1;

int oldX=0, oldY=0;
float rX=15, rY=0;
int state =1 ;
float dist=-23;
const int GRID_SIZE=10;
float pointSize = 10;
float spacing =  float(sizeX)/(numX+1);                         // Spacing of particles

glm::vec3 gravity=glm::vec3(0.0f,-9.81f,0.0f);

GLint viewport[4];
GLdouble MV[16];
GLdouble P[16];

glm::vec3 Up=glm::vec3(0,1,0), Right, viewDir;

struct timeval t1, t2;           // ticks
double frameTimeQP=0;
float frameTime =0 ;

float startTime =0, fps=0;
int totalFrames=0;

int whichIndex = 0;
char info[256]={0};

glm::mat3       I=glm::mat3(1);         //identity matrix
float nu =  0.4f;               //Poisson ratio
float Y = 300.0f;            //Young modulus
float density = 10000.f;        //material density
float kv=100, damping=50.0f;   //constant used in Eq. 22 page 5
float d15 = Y / (1.0f + nu) / (1.0f - 2 * nu);
float d16 = (1.0f - nu) * d15;
float d17 = nu * d15;
float d18 = Y / 2 / (1.0f + nu);
glm::vec3 D(d16, d17, d18);

float scFac = 0; //scaling constant

vector<Phyxel*> phyxels;

void OnRender();
void OnReshape(int nw, int nh);
void OnIdle();
void OnMouseMove(int x, int y);
void OnMouseDown(int button, int s, int x, int y);
void OnKey(unsigned char k,int , int);
void InitGL();

void GetKNearestNeighbors(int index, int k, vector<float>& dis, vector<Neighbor>& n);
void ComputeRadiusAndSupportRadii(vector<float> &dists, float& r, float& h);
void FillNeighs(int index);
void ComputeScalingConstant();
void ComputeMass(float dm, int index);

void ComputeDensityAndVolume(int index);
void ComputeInvMomentMatrix(int index);

void StepPhysics(float dt);
void UpdateForces();
void ComputeJacobians();
void ComputeStrainAndStress();
void IntegrateLeapFrog(float dt);

int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
    glutInitWindowSize(width, height);
    glutCreateWindow("DDD ---- ");
    
    glutDisplayFunc(OnRender);
    glutReshapeFunc(OnReshape);
    glutIdleFunc(OnIdle);
    
    glutMouseFunc(OnMouseDown);
    glutMotionFunc(OnMouseMove);
    
    glutKeyboardFunc(OnKey);
    glewInit();
    
    // Initialize
    InitGL();
    
    glutMainLoop();
    
    return 0;
}

void OnRender() {
    int i=0;
    float newTime = (float) glutGet(GLUT_ELAPSED_TIME);
    frameTime = newTime-currentTime;
    currentTime = newTime;
    
    gettimeofday(&t2, NULL);
    frameTimeQP = (t2.tv_sec - t1.tv_sec) * 1000000.0;
    t1=t2;
    accumulator += frameTimeQP;
    
    ++totalFrames;
    
    if((newTime-startTime)>100)
    {
        float elapsedTime = (newTime-startTime);
        fps = (totalFrames/ elapsedTime)*100000 ;
        startTime = newTime;
        totalFrames = 0;
    }
    
    glutSetWindowTitle("DDD");
    glClear(GL_COLOR_BUFFER_BIT| GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    glTranslatef(0,0,dist);
    glRotatef(rX,1,0,0);
    glRotatef(rY,0,1,0);
    
    glGetDoublev(GL_MODELVIEW_MATRIX, MV);
    viewDir.x = (float)-MV[2];
    viewDir.y = (float)-MV[6];
    viewDir.z = (float)-MV[10];
    Right = glm::cross(viewDir, Up);
    
    
    //draw points
    glBegin(GL_POINTS);
    for(i=0;i<num_particle;i++) {
        glm::vec3 p = phyxels[i]->X;
        /*
        if(i==selected_index)
            glColor3f(0,1,1);
        else
            glColor3f(1,0,0);
            */
        glColor3f(1,1,1);
        glVertex3f(p.x,p.y,p.z);
    }
    
    glEnd();
    
    glutSwapBuffers();
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
    
    if ( accumulator >= timeStep )
    {
        StepPhysics(timeStep);
        accumulator -= timeStep;
    }
    
    glutPostRedisplay();
}

void OnMouseDown(int button, int s, int x, int y)
{
    if (s == GLUT_DOWN)
    {
        oldX = x;
        oldY = y;
        int window_y = (height - y);
        float norm_y = float(window_y)/float(height/2.0);
        int window_x = x ;
        float norm_x = float(window_x)/float(width/2.0);
        
        float winZ=0;
        glReadPixels(x, height-y, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ);
        if(winZ==1)
            winZ=0;
        double objX=0, objY=0, objZ=0;
        gluUnProject(window_x,window_y, winZ, MV, P, viewport, &objX, &objY, &objZ);
        glm::vec3 pt(objX, objY, objZ);
        int i=0;
        for(i=0;i<num_particle;i++) {
            if( glm::distance(phyxels[i]->X,pt)<1) {
                selected_index = i;
                printf("Intersected at %d\n",i);
                break;
            }
        }
    }
    
    if(button == GLUT_MIDDLE_BUTTON)
        state = 0;
    else
        state = 1;
    
    if(s==GLUT_UP) {
        selected_index= -1;
        glutSetCursor(GLUT_CURSOR_INHERIT);
    }
}

void OnMouseMove(int x, int y)
{
    if(selected_index == -1) {
        if (state == 0)
            dist *= (1 + (y - oldY)/60.0f);
        else
        {
            rY += (x - oldX)/5.0f;
            rX += (y - oldY)/5.0f;
        }
    } else {
        float delta = 1500/abs(dist);
        float valX = (x - oldX)/delta;
        float valY = (oldY - y)/delta;
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
    oldX = x;
    oldY = y;
    
    glutPostRedisplay();
}

void OnKey(unsigned char k,int , int) {
    switch(k) {
        case 'a':whichIndex--;break;
        case 'd':whichIndex++;break;
        case ',':Y-=500;break;
        case '.':Y+=500;break;
        case '[':nu-=0.01f;break;
        case ']':nu+=0.01f;break;
    }

    if(nu>0.49999f)     nu=0.49f;
    if(nu<0)            nu=0;
    if(Y<0.01f)         Y=0.01f;
    if(Y>175000000)     Y=175000000;
    if(whichIndex<0)    whichIndex = 0;
    
    /* Update D */
    d15 = Y / (1.0f + nu) / (1.0f - 2 * nu);
    d16 = (1.0f - nu) * d15;
    d17 = nu * d15;
    d18 = Y / 2 / (1.0f + nu);

    D.x=d16;
    D.y=d17;
    D.z=d18;

    whichIndex = (whichIndex%num_particle);

    glutPostRedisplay();
    
}

void InitGL() {
    startTime = (float)glutGet(GLUT_ELAPSED_TIME);
    currentTime = startTime;
    
    gettimeofday(&t1, NULL);
    
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_POINT_SMOOTH);
    int i = 0, j = 0, k = 0, count=0;
    
    float ypos = 4.0f;

    for (i = 0; i < num_particle; i++) {
        Phyxel* phyxel = new Phyxel();
        phyxels.push_back(phyxel);
    }
    
    for(k = 0; k < numZ; k++) {
        for(j = 0; j < numY; j++) {
            for(i = 0; i < numX; i++) {
                phyxels[count++]->setX(glm::vec3(((float(i)/(numX-1)) ) *  sizeX,
                                       ((float(j)/(numY-1))*2-1)* hsizeY + ypos,
                                       ((float(k)/(numZ-1))*2-1)* hsizeZ));
                phyxels[count-1]->setXi(phyxels[count-1]->X);
                if (phyxels[count-1]->X.y == ypos + hsizeY){
                    phyxels[count-1]->isFixed = true;
                }
                
            }
        }
    }
    
    for(i = 0; i < num_particle; ++i) {
        vector<float> dist;
        GetKNearestNeighbors(i, 10, dist, phyxels[i]->neighbors);
        ComputeRadiusAndSupportRadii(dist, phyxels[i]->r, phyxels[i]->h);
    }
    
    for(i = 0; i < num_particle; i++) {
        FillNeighs(i);
    }
    
    ComputeScalingConstant();
    
    for(i=0; i < num_particle; ++i)
        ComputeMass(density, i);
    
    for(i=0; i < num_particle; ++i)
        ComputeDensityAndVolume(i);
    
    for(i=0; i < num_particle; ++i) {
        ComputeInvMomentMatrix(i);
    }

    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    //glPolygonMode(GL_BACK, GL_LINE);
    glPointSize(pointSize );
    
    //wglSwapIntervalEXT(0);
}


typedef std::pair<int, float> mypair;

struct Cmp {
    bool operator()(const mypair &lhs, const mypair &rhs) {
        return lhs.second < rhs.second;
    }
};

void GetKNearestNeighbors(int index, int k, vector<float>& dis, vector<Neighbor>& n) {
    
    vector<mypair> distances;

    //Get the distances of current point to all other points
    for(int i = 0; i < num_particle; i++) {
        if(index!=i) {
            mypair m;
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

void ComputeRadiusAndSupportRadii(vector<float> &dists, float& r, float& h)
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

void FillNeighs(int index)
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

void ComputeScalingConstant() {
    scFac = 0.f;
    int i=0;
    for ( i =0; i < num_particle; i++ ) {
        float sum = 0.f;

        vector<Neighbor>& pNeighbor = phyxels[i]->neighbors;
        for(size_t j = 0; j < pNeighbor.size(); j++) {
            sum += pow(phyxels[pNeighbor[j].j]->r, 3) * pNeighbor[j].w;
        }
        scFac += 1.0f / sum;
    }
    // This is the common scaling factor to compute the mass of each phyxel.
    // See last paragraph of Section 3.2 on page 4
    scFac /= num_particle;

    printf("Scaling factor: %3.3f\n", scFac);
}

void ComputeMass(float dm, int index)
{
    // See last paragraph of Section 3.2 on page 4
    phyxels[index]->m = scFac * pow(phyxels[index]->r, 3) * dm;
}

void ComputeDensityAndVolume(int index)
{
    // See last paragraph of Section 3.2 on page 4
    phyxels[index]->rho = 0.f;
    vector<Neighbor>& pNeighbor = phyxels[index]->neighbors;
    for(size_t i = 0; i < pNeighbor.size(); i++)
        phyxels[index]->rho += phyxels[pNeighbor[i].j]->m * pNeighbor[i].w;       // Eq. 10 page 4
    phyxels[index]->vol =  phyxels[index]->m / phyxels[index]->rho;
}

void ComputeInvMomentMatrix(int index)
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

void StepPhysics(float dt) {
    UpdateForces();
    IntegrateLeapFrog(dt);
}

void UpdateForces()
{
    //Calculate external force
    for(int i=0;i<num_particle;i++) {
        if(!phyxels[i]->isFixed)
            phyxels[i]->F = glm::vec3(0,gravity.y,0);
        else
            phyxels[i]->F = glm::vec3(0);

        //Add velocity damping
        phyxels[i]->F -= phyxels[i]->v * damping;
    }

    ComputeJacobians();

    ComputeStrainAndStress();


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

void ComputeJacobians()
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

void ComputeStrainAndStress()
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

void IntegrateLeapFrog(float dt) {
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
    UpdateForces();

    //V_(i+1) = V_i + ((a_i+a_(i+1))/2)*dt^2
    for(i=0;i<num_particle;i++) {
        if(!phyxels[i]->isFixed)
            phyxels[i]->v += ((phyxels[i]->F/phyxels[i]->m 
                + phyxels[i]->acc)*half_dt2);
    }

}