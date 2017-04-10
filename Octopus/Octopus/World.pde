// A world consist of a spatial hash map
import java.util.Map;
import java.util.List;
import java.util.Collections;

class Vector3D
{
  int x, y, z;
  public int distance2(Vector3D aVector)
  {
    return int (pow((aVector.x-x),2) + pow((aVector.y-y),2) + pow((aVector.z-z),2));
  }
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
  
  // Set the bin size to the radius size!!!!!!
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
  }
  
  private int hash(Vector3D position)
  {
    //return (position.x * p1) ^ (position.y * p2) ^ (position.z * p3) % bucketNum;
    return position.x + position.y * worldSize.x + position.z * worldSize.x * worldSize.y;
  }
  
  public boolean add(Phyxel p)
  {
    hashMap.get(hash(p.matCoord)).add(p);
    return true;
  }
  
  public ArrayList<Phyxel> get(int cell)
  {
    return hashMap.get(cell);
  }
  

  
  public ArrayList<Phyxel>getNearest(Vector3D position, int num)
  {
    ArrayList<Phyxel> result = new ArrayList<Phyxel>();
    ArrayList<Phyxel> cellElements = hashMap.get(hash(position));

    int iterNum = cellElements.size();
    int cellSize2 = cellSize * cellSize;
    for(int i = 0; i < iterNum; i++)
    {
      Phyxel currentElement = cellElements.get(i);
      int distance = position.distance2(currentElement.position);
      // remove itself
      if(distance < cellSize2 && distance != 0)
      {
        result.add(currentElement);
      }
    }
    
    // TODO: sort and get the smallest num elements

    
    return result;
  }

}