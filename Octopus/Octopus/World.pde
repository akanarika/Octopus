// A world consist of a spatial hash map
import java.util.Map;
import java.util.List;

class Vector3D
{
  float x, y, z;
}

class World
{
  static final int p1 = 73856093;
  static final int p2 = 19349663;
  static final int p3 = 83492791;
  
  HashMap<Integer, ArrayList<Phyxel>> hashMap;
  public int cellSize;
  public Vector3D worldSize;
  public int bucketNum;
  public Vector3D center;
  
  public World(Vector3D _center, Vector3D _worldSize, int _cellSize, int _bucketNum)
  {
    center = _center;
    worldSize = _worldSize;
    cellSize = _cellSize;
    bucketNum = _bucketNum;
    
    hashMap = new HashMap <Integer, ArrayList<Phyxel>>(bucketNum);
    for(int i=0; i<bucketNum; i++)
    {
      hashMap.put(i, new ArrayList<Phyxel>());
    }
  
  private int hash(Vector3D position)
  {
    return (position.x * p1) ^ (position.y * p2) ^ (position.z * p3) % bucketNum;
  }
  
  public bool add(Phyxel p)
  {
    hashMap.get(hash(p.matCoord)).add(p);
    return true;
  }
  
  public ArrayList<Phyxel> get(int cell)
  {
    return hashMap.get(cell);
  }

}