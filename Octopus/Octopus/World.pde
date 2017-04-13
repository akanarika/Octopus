// A world consist of a spatial hash map
import java.util.Map;
import java.util.List;
import java.util.Collections;

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
  public int distance2(Vector3D aVector)
  {
    return int (pow((aVector.x-x),2) + pow((aVector.y-y),2) + pow((aVector.z-z),2));
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
  
  ArrayList<ArrayList<Float>> wlist;
  ArrayList<ArrayList<Vector3D>> xlist;
   
  public Vector3D worldSize;
  private Vector3DI worldSizeInCell;
  private float cellSize;
  private float searchRadius;
  public Vector3D origin;
  public Phyxel[] phyxels;
  
  /* w_ij = W(|x_j - x_i|, h_i)
   * W(r, h) is the polynomial kernel
   */
  ArrayList<ArrayList<Float>> w;
  /* x_ij = x_j - x_i
  */
  ArrayList<ArrayList<Vector3D>> x;
    
  public World(Vector3D _origin, Vector3D _worldSize, float _searchRadius, PShape _model)
  {
    origin = _origin;
    worldSize = _worldSize;
    searchRadius = _searchRadius;
    cellSize = 2 * _searchRadius; // Set the bin size to the 2 * radius size!!!!!!
    
    worldSizeInCell = new Vector3DI (ceil(worldSize.x / cellSize), ceil(worldSize.y / cellSize), ceil(worldSize.z / cellSize));
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
  
  // Calculate h of each phyxel
  public void initH()
  {
    for (int i = 0; i < phyxels.length; i++)
    {
      ArrayList<Integer> neighbours = phyxels[i].getNeighbours();
      float average_r = 0;
      for (int ii = 0; ii < neighbours.size(); ii++)
      {
        average_r += sqrt((phyxels[i].matCoord).distance2(phyxels[ii].matCoord));
      }
      average_r /= neighbours.size();
      float h = average_r * 3;
      phyxels[i].h = h;
    }
  }
  
  // Fill w_ij list and x_ij list and density and volume
  public void initXWDV() 
  { 
    for (int i = 0; i<phyxels.length; i++) 
    {
      ArrayList<Float> _wlist = new ArrayList<Float>();
      ArrayList<Vector3D> _xlist = new ArrayList<Vector3D>();
      float rho = 0;
      float m_j = phyxels[i].mass;
      for (int j =0; j < phyxels.length; j++)
      {
        float w = 0;
        
        Vector3D x_ij = (phyxels[j].matCoord).subtract(phyxels[i].matCoord); //  x_ij = x_j - x_i
        _xlist.add(x_ij);
        float r = x_ij.magnitude();  // r = |x_j - x_i|
        float h = phyxels[i].h;
        
        if (r >= h) w = 0;
        else w = 315 * pow((h*h - r*r), 3) / (64 * (float)Math.PI * pow(h, 9));
        _wlist.add(w);
        rho += m_j * w;
      }
      phyxels[i].density = rho;
      phyxels[i].volume = phyxels[i].mass / rho;
      wlist.add(_wlist);
      xlist.add(_xlist);
    }
  }

}