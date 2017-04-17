class Vector3D
{
  float x, y, z;
  
  public Vector3D(float _x, float _y, float _z)
  {
    x = _x;
    y = _y;
    z = _z;
  }
  
  /* Return the power of the distance between two vectors*/
  public float distance2(Vector3D aVector)
  {
    return  (pow((aVector.x-x),2) + pow((aVector.y-y),2) + pow((aVector.z-z),2));
  }
  
  /* Return the distance between two vectors*/
  public float distance(Vector3D aVector)
  {
    return sqrt(pow((aVector.x-x),2) + pow((aVector.y-y),2) + pow((aVector.z-z),2));
  }
  
  public float magnitude()
  {
    return sqrt(x*x + y*y + z*z);
  }
  
  public Vector3D subtract(Vector3D aVector)
  {
    return new Vector3D(x - aVector.x, y - aVector.y, z - aVector.z);
  }
}


class Vector3DI
{
  int x, y, z;
  public Vector3DI(int _x, int _y, int _z){ x = _x; y = _y; z = _z;}
  public int total() { return x * y * z; }
}

class IdxDist2Pair
{
  int index;
  float dist2;
  public IdxDist2Pair(int _i, float _dist2)
  {
    index = _i;
    dist2 = _dist2;
  }
}