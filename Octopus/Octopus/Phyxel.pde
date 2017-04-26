final float materialDensity = 1.0;  // should change
final float densityScalar = 1.0;  // should change
final int NUM_OF_NEAREST_NEIGHBOUR = 10;

class Phyxel
{
  ArrayList<Integer> neighbours;
  float mass;
  float density;      // rho = sum_j(m_j * w_ij)
  float volume;       // v_i = m_i / rho_i
  float h;
  Vector3D matCoord;  // x
  Vector3D v;
  //Vector3D position;
  Vector3D u;
  Vector3D f;
  int index;          // The index in the world array    
  
  public Phyxel()
  {
    // init mass
    // init hi
      // set nearest neighbours
      // compute average r
      // hi = 3(>1) * r
  }
  
  public Phyxel(Vector3D _matCoord, int i)
  {
    matCoord = _matCoord;
    index = i;
    u = new Vector3D(0,0,0);
    v = new Vector3D(0,0,0);
  }
  
  public void setNeighbors(ArrayList<Integer> _neighbours) 
  {
    neighbours = new ArrayList<Integer>(_neighbours);
  }
  
  public ArrayList<Integer> getNeighbours()
  {
    return neighbours;
  }
  
  
  void update()
  {
    // calculate A
    // calculate derivatives
    // jacobian, stress, strain
    // force
    // acceleration
    // integrate
  }
   
}