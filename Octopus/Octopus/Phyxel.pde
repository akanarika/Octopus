final float materialDensity = 1.0;  // should change
final float densityScalar = 1.0;  // should change
final int NUM_OF_NEAREST_NEIGHBOUR = 10;

class Phyxel
{
  Phyxel[] neighbours;
  float mass;
  float density;
  Vector3D matCoord;
  Vector3D position;
  
  public Phyxel()
  {
    // init mass
    // init hi
      // set nearest neighbours
      // compute average r
      // hi = 3(>1) * r
  }
  
  void findNearestNeighbours()
  {
    //Spatial Hashing
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