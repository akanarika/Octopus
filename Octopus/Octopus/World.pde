// A world consist of a spatial hash map
import java.util.Map;
import java.util.List;
import java.util.Collections;

class Vector3D
{
  float x, y, z;
  
  /* Return the power of the distance between two vectors*/
  public int distance2(Vector3D aVector)
  {
    return int (pow((aVector.x-x),2) + pow((aVector.y-y),2) + pow((aVector.z-z),2));
  } 
}

class Vector3DI
{
  int x, y, z;
  public Vector3DI(int _x, int _y, int _z){ x = _x; y = _y; z = _z;}
  public int total() { return x * y * z; }
}

/* World : seperate the space into small squares
 * y
 * ^    z
 * |   /
 * |  /
 * | /
 * | - - - - -> x
 * 0
 */
class World
{
  HashMap<Integer, ArrayList<Integer>> hashMap;
   
  public Vector3D worldSize;
  private Vector3DI worldSizeInCell;
  private float cellSize;
  private float searchRadius;
  public Vector3D origin;
  public Phyxel[] phyxels;
    
  public World(Vector3D _origin, Vector3D _worldSize, float _searchRadius, PShape _model)
  {
    origin = _origin;
    worldSize = _worldSize;
    searchRadius = _searchRadius;
    cellSize = 2 * _searchRadius; // Set the bin size to the 2 * radius size!!!!!!
    
    worldSizeInCell = new Vector3DI (ceil(worldSize.x / cellSize), ceil(worldSize.y / cellSize), ceil(worldSize.z / cellSize);
    int bucketNum = worldSizeInCell.total();
    
    hashMap = new HashMap <Integer, ArrayList<Integer>>(bucketNum);
    for(int i=0; i<bucketNum; i++)
    {
      hashMap.put(i, new ArrayList<Integer>());
    }
    
    // Load model into Phyxels
    initModel(_model);
  }
  
  void initModel(PShape model)
  {
   
  }
  
  public int toIndex(Vector3D position)
  {
    Vector3DI cellIndex=  coordToCellIndex(position);
    return cellIndexToIndex(cellIndex);
  }
    // cell: [Min, Max)
  private Vector3DI coordToCellIndex(Vector3D worldPos)
  {
    return new Vector3DI(floor(worldPos.x / worldSizeInCell.x) , floor(worldPos.y / worldSizeInCell.y), floor(worldPos.z/worldSizeInCell.z));
  }
  
  /* hash a cell_space index into an integer */
  private int cellIndexToIndex(Vector3DI pos)
  {
    //return (position.x * p1) ^ (position.y * p2) ^ (position.z * p3) % bucketNum;
    return pos.x + pos.y * worldSizeInCell.x + pos.z * worldSizeInCell.y * worldSizeInCell.x;
  }
  
  public boolean add(Phyxel p)
  {
    hashMap.get( toIndex(p.position) ).add(p.index);
    return true;
  }
  
  public ArrayList<Integer> getListAt(int cell)
  {
    return hashMap.get(cell);
  }
  
  public void clear(int cell)
  {
    hashMap.put (cell, new ArrayList<Integer>());
  }
  
  public ArrayList<Phyxel>getNearest(Vector3D position, int num)
  {
    ArrayList<Phyxel> result = new ArrayList<Phyxel>();
   
     
    ArrayList<Phyxel> cellElements = hashMap.get(hash(position));

    int iterNum = cellElements.size();
    float cellSize2 = cellSize * cellSize;
    
    for(int i = 0; i < iterNum; i++)
    {
      Phyxel currentElement = cellElements.get(i);
      int distance2 = position.distance2(currentElement.position);
      // remove itself
      if(distance2 < cellSize2 && distance2 != 0)
      {
        result.add(currentElement);
      }
    }
    
    // TODO: sort and get the smallest num elements

    
    return result;
  }

}