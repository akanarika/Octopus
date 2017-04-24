// A world consist of a spatial hash map
import java.util.Map;
import java.util.List;
import java.util.Collections;
import java.util.Comparator;
import java.util.PriorityQueue;
import Jama.*;

  /* World : seperate the space into small squares
   * y
   * ^    z
   * |   /
   * |  /
   * | /
   * | - - - - -> x
   * 0
   */
  class World {
    // Spatial Hashing Map
    HashMap<Integer, ArrayList<Integer>> hashMap;
    public Phyxel[] phyxels;

    // World variables
    Vector3D worldSize;
    Vector3DI worldSizeInCell;
    float cellSize;
    float searchR;                // search radius
    
    // Calculation Constants
    ArrayList<ArrayList<Float>> wlist;  // w_ij = W(|x_j - x_i|, h_i), W(r, h) is the polynomial kernel
    ArrayList<ArrayList<Vector3D>> xlist;  // x_ij = x_j - x_i
    ArrayList<WB_M33> Alist;  // A_i = sum_j(x_ij * x_ij^T * w_ij) //<>//

        
    // Using HE_Mesh because it has no repeat vertex
    // Force the origin to be (0, 0, 0) and put the model in the positive area

    public World(Vector3D _worldSize, float _searchRadius, HE_Mesh _model)
    { 
    
      worldSize = _worldSize;
      searchR = _searchRadius;
      cellSize = 2 * _searchRadius;   // Setting the bin size to the 2 * radius size
      worldSizeInCell = new Vector3DI (ceil(worldSize.x / cellSize), ceil(worldSize.y / cellSize), ceil(worldSize.z / cellSize));
      
      int bucketNum = worldSizeInCell.total();
      hashMap = new HashMap <Integer, ArrayList<Integer>>(bucketNum);
      for(int i=0; i<bucketNum; i++){
        hashMap.put(i, new ArrayList<Integer>());
      }
      
      // Load model into Phyxels
      _loadModel(_model);
      
    }
    
    // Load the model into phyxels. 
    void _loadModel(HE_Mesh model)
    {
      WB_Coord [] vertices = model.getVerticesAsArray();
      //List<WB_Coord> vertexes = model.getPoints();
      int vertCount = vertices.length;
      print(vertCount);
      phyxels = new Phyxel[vertCount];
      for (int i = 0; i < vertCount; i ++)
      {
        // Load into Phyxel
        WB_Coord vertPos = vertices[i];
        Phyxel newPhyxel = new Phyxel( new Vector3D(vertPos.xf(), vertPos.yf(), vertPos.zf()) , i );
        phyxels[i] =  newPhyxel;
        
        // Add Phyxel to the hash map
        AddToCells(newPhyxel);
      }
    }
    
    // Convert a 3D matCoord to an array index
    public int CoordToIndex(Vector3D matCoord)
    {
      // Check boundary
      if(matCoord.x < 0 || matCoord.y < 0 ||matCoord.z < 0)
        return -1;
        
      Vector3DI cellIndex = _coordToCellIndex(matCoord);
      return _cellIndexToIndex(cellIndex);
    }
    
    // cell: [Min, Max)
    Vector3DI _coordToCellIndex(Vector3D worldPos)
    {
      return new Vector3DI(floor(worldPos.x / worldSizeInCell.x) , floor(worldPos.y / worldSizeInCell.y), floor(worldPos.z/worldSizeInCell.z));
    }
    
    // hash a cell_space index into an integer
    int _cellIndexToIndex(Vector3DI pos)
    {
      return pos.x + pos.y * worldSizeInCell.x + pos.z * worldSizeInCell.y * worldSizeInCell.x;
    }
    
    public void AddToCells(Phyxel p)
    {
      ArrayList<Integer> cellIds = GetIdsForPhyxel(p);
      for(int i=0; i<cellIds.size(); i++)
      {
        ArrayList<Integer> cell = hashMap.get(cellIds.get(i));
        AddToList(cell, p.index);
      }
    }
    
    // get Phyxel p's neighbourhood bucket 
    public ArrayList<Integer>GetIdsForPhyxel(Phyxel p)
    {
      ArrayList<Integer> cellIdsToAdd = new ArrayList<Integer>();
      Vector3D matCoord = p.matCoord;
      
      int index;
      // Bottom
      index = CoordToIndex(new Vector3D(matCoord.x, matCoord.y - searchR, matCoord.z));
      if(index != -1){
        AddToList(cellIdsToAdd, index);
      }
      // Top
      index = CoordToIndex(new Vector3D(matCoord.x, matCoord.y + searchR, matCoord.z));
      if(index != -1){
        AddToList(cellIdsToAdd, index);
      }
      //Left
      index = CoordToIndex(new Vector3D(matCoord.x - searchR, matCoord.y, matCoord.z));
      if(index != -1)
      {
        AddToList(cellIdsToAdd, index);
      }
      // Right
      index = CoordToIndex(new Vector3D(matCoord.x + searchR, matCoord.y, matCoord.z));
      if(index != -1)
      {
        AddToList(cellIdsToAdd, index);
      }
      //Front
      index = CoordToIndex(new Vector3D(matCoord.x, matCoord.y, matCoord.z - searchR));
      if(index != -1)
      {
        AddToList(cellIdsToAdd, index);
      }
      // Back
      index = CoordToIndex(new Vector3D(matCoord.x, matCoord.y, matCoord.z + searchR));
      if(index != -1)
      {
        AddToList(cellIdsToAdd, index);
      }
      
      return cellIdsToAdd;
    }
    
    // Wrapper to add unrepeated element to the list
    public void AddToList(ArrayList<Integer>list, int element)
    {
      if(!list.contains((element))){
        list.add(element);
      }
    }
    
    public ArrayList<Integer> GetListAt(int cell)
    {
      return hashMap.get(cell);
    }
    
    //public void ClearAt(int cell)
    //{
    //  hashMap.put (cell, new ArrayList<Integer>());
    //}
    
    public void ClearAll()
    {
      hashMap.clear();
      int total = worldSizeInCell.total();
      for(int i=0; i<total; i++)
      {
        hashMap.put (i, new ArrayList<Integer>());
      }
    }
    
    public Phyxel FindPhyxelAtLocation(Vector3D pos)
    {
      int index = CoordToIndex(pos);
      ArrayList<Integer> nearbyVerts = GetListAt(index);
      int iterNum = nearbyVerts.size();
      for(int i=0; i < iterNum; i++)
      {
        Vector3D vert = phyxels[i].matCoord;
        if(vert.x == pos.x && vert.y == pos.y && vert.z == pos.z)
          return phyxels[i];
      }
      return null;
    }
    
    // Get num neighbours for pixel p
    public ArrayList<Integer>GetNeighbours(Phyxel p, int num)
    {
     
      ArrayList<Integer> result = new ArrayList<Integer>();
      ArrayList<Integer> nearbyVerts = new ArrayList<Integer>();
      ArrayList<Integer> cellIds = GetIdsForPhyxel(p);
      Vector3D vertPos = p.matCoord;
      
      // Get Nearby Vertices
      int iterNum;
      iterNum = cellIds.size();
      for(int i=0; i<iterNum; i++)
      {
        int cellId = cellIds.get(i);
        nearbyVerts.addAll(GetListAt(cellId));
      }
      
      // minHeap
      PriorityQueue<IdxDist2Pair> pq = new PriorityQueue<IdxDist2Pair>(new Comparator<IdxDist2Pair>()
      {
        public int compare(IdxDist2Pair p1, IdxDist2Pair p2){
          if(p1.dist2 < p2.dist2) return -1;
          if(p1.dist2 > p2.dist2) return 1;
          return 0;
        }
      });
      
      // Iterate Over Every Point
      float searchR2 = searchR * searchR;
      iterNum = nearbyVerts.size();
      for(int i=0; i<iterNum; i++)
      {
        int nearbyIndex = nearbyVerts.get(i);
        Vector3D nearbyPos = phyxels[nearbyIndex].matCoord;
        float distance2 = vertPos.distance2(nearbyPos);
        //print(nearbyIndex + " " + nearbyPos.x + " " + nearbyPos.y + " " + nearbyPos.z + " " + distance2 + "\n");
        
        if(distance2 < searchR2){
          if(nearbyIndex != p.index){
            IdxDist2Pair pair = new IdxDist2Pair(nearbyIndex, distance2);
            pq.offer(pair);
          }
        }
      }
      
      iterNum = min(num, pq.size());
      int count = 0;
      while(count < iterNum)
      {
        
        IdxDist2Pair pair = pq.poll();
        if (pair==null) break;
        if(!result.contains(pair.index)){
          result.add(pair.index);
          count ++;
        }
        
      }
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
    
    // Calculate matrix A
    // A_i = sum_j(x_ij * x_ij^T * w_ij)
    public void initA()
    {
      WB_M33 A = new WB_M33(0,0,0,0,0,0,0,0,0);
      WB_M33 a;
      for(int i=0; i < phyxels.length; i++)
      {
        ArrayList<Integer> neighbours = phyxels[i].getNeighbours();
        for(int j=0; j < neighbours.size(); j++)
        {
          a = new WB_M33(xlist.get(i).get(j).x * xlist.get(i).get(j).x, 
                         xlist.get(i).get(j).x * xlist.get(i).get(j).y,
                         xlist.get(i).get(j).x * xlist.get(i).get(j).z,
                         xlist.get(i).get(j).y * xlist.get(i).get(j).x, 
                         xlist.get(i).get(j).y * xlist.get(i).get(j).y,
                         xlist.get(i).get(j).y * xlist.get(i).get(j).z,
                         xlist.get(i).get(j).z * xlist.get(i).get(j).x, 
                         xlist.get(i).get(j).z * xlist.get(i).get(j).y,
                         xlist.get(i).get(j).z * xlist.get(i).get(j).z);
          A.add(a);
          A.mul(wlist.get(i).get(j));
        }
        Alist.add(A);
      }
    }
    
    //  Force of each phyxel has to be initialized to 0 at each frame
    //  (which has not been done yet!!!)
    public void updateF()
    {
      Vector3D _u, _v, _w;
      Vector3D du, dv, dw;
      WB_M33 U, UT, UUT, J;  // deriv U; deriv U transpose; Jacob
      WB_M33 epsilon = new WB_M33();
      WB_M33 sigma = new WB_M33();  // ε, σ
      WB_M33 Fe = new WB_M33();
      WB_M33 Fv = new WB_M33();
      WB_M33 Fe_add_Fv = new WB_M33();
      WB_Vector Ju, Jv, Jw;
      WB_Vector JvJw, JwJu, JuJv;
      WB_Vector di, dj;
      WB_Vector _fi = new WB_Vector();
      WB_Vector _fj = new WB_Vector();
      Vector3D _di = new Vector3D(0,0,0);  // -sum_j(x_ij * w_ij)
      float C = 4;  // C
      float k = 1; // kv ??????
      for (int i=0; i < phyxels.length; i++)
      {
        _u = new Vector3D(0,0,0);
        _v = new Vector3D(0,0,0);
        _w = new Vector3D(0,0,0);
        ArrayList<Integer> neighbours = phyxels[i].getNeighbours();
        for (int j=0; j < neighbours.size(); j++)
        {
          _u.add(xlist.get(i).get(j).mult((phyxels[j].u.x - phyxels[i].u.x) * wlist.get(i).get(j)));
          _v.add(xlist.get(i).get(j).mult((phyxels[j].u.y - phyxels[i].u.y) * wlist.get(i).get(j)));
          _w.add(xlist.get(i).get(j).mult((phyxels[j].u.z - phyxels[i].u.z) * wlist.get(i).get(j)));
          _di.add((xlist.get(i).get(j)).mult(wlist.get(i).get(j)));
        }
        du = new Vector3D(Alist.get(i).m11 * _u.x + Alist.get(i).m12 * _u.y + Alist.get(i).m13 * _u.z,        
                          Alist.get(i).m21 * _u.x + Alist.get(i).m22 * _u.y + Alist.get(i).m23 * _u.z,
                          Alist.get(i).m31 * _u.x + Alist.get(i).m32 * _u.y + Alist.get(i).m33 * _u.z);
        dv = new Vector3D(Alist.get(i).m11 * _v.x + Alist.get(i).m12 * _v.y + Alist.get(i).m13 * _v.z,        
                          Alist.get(i).m21 * _v.x + Alist.get(i).m22 * _v.y + Alist.get(i).m23 * _v.z,
                          Alist.get(i).m31 * _v.x + Alist.get(i).m32 * _v.y + Alist.get(i).m33 * _v.z);
        dw = new Vector3D(Alist.get(i).m11 * _w.x + Alist.get(i).m12 * _w.y + Alist.get(i).m13 * _w.z,        
                          Alist.get(i).m21 * _w.x + Alist.get(i).m22 * _w.y + Alist.get(i).m23 * _w.z,
                          Alist.get(i).m31 * _w.x + Alist.get(i).m32 * _w.y + Alist.get(i).m33 * _w.z);
        // Deriv u
        U = new WB_M33(du.x, dv.x, dw.x, du.y, dv.y, dw.y, du.z, dv.z, dw.z);
        
        // u transpose
        UT = U;
        UT.transpose();
        
        // Jacobi
        J = new WB_M33(du.x + 1, du.y, du.z, dv.x, dv.y + 1, dv.z, dw.x, dw.y, dw.z + 1);
        
        // U * UT
        UUT = WB_M33.mul(U, UT);
        
        // U = U + UT
        U.add(UT);
        // epsilon = U + UT + U * UT
        U.addInto(UUT, epsilon);
        // sigma = C * epsilon
        sigma = epsilon;
        sigma.mul(C);
        
        Ju = J.row(0);
        Jv = J.row(1);
        Jw = J.row(2);
        JvJw = Jv.cross(Jw);
        JwJu = Jv.cross(Jw);
        JuJv = Ju.cross(Jv);
        di = WB_M33.mulToVector(Alist.get(i).inverse(), _di.to_WB_Vector());  //  di = A^-1*(-sum_j(x_ij*w_ij))
        WB_M33.mulInto(J, sigma, Fe);
        Fe.mul(-2 * phyxels[i].volume);
        Fv = new WB_M33(JvJw.coords()[0], JvJw.coords()[1], JvJw.coords()[2],
                        JwJu.coords()[0], JwJu.coords()[1], JwJu.coords()[2],
                        JuJv.coords()[0], JuJv.coords()[1], JuJv.coords()[2]);
        Fv.mul(-phyxels[i].volume * k * (J.det() - 1));
        Fe.addInto(Fv, Fe_add_Fv);
        _fi = WB_M33.mulToVector(Fe_add_Fv, di);
        phyxels[i].f.add(new Vector3D(_fi.coords()));
        for (int j=0; j < neighbours.size(); j++)
        {
          dj = WB_M33.mulToVector(Alist.get(i).inverse(), (xlist.get(i).get(j)).mult(wlist.get(i).get(j)).to_WB_Vector());
          _fj = WB_M33.mulToVector(Fe_add_Fv, dj);
          phyxels[j].f.add(new Vector3D(_fj.coords()));
        }
      }
      
    }
    
    // Display the nearest neighbours on the static model
    public void draw(int idx)
    {
      int iterNum = phyxels.length;
      stroke(0, 0, 0);
      for(int i=0; i < iterNum; i++)
      {
        Vector3D matCoord = phyxels[i].matCoord;
        vertex(matCoord.x, matCoord.y, matCoord.z);
      }
      
      strokeWeight(10);
      stroke(0,255,0);
      vertex(phyxels[idx].matCoord.x, phyxels[idx].matCoord.y, phyxels[idx].matCoord.z);
      ArrayList< Integer > neighbour = GetNeighbours(phyxels[idx],10);
      iterNum = neighbour.size();
      //print(iterNum);
      
      
      strokeWeight(5);
      stroke(255, 0, 0);
      for(int i=0; i < iterNum; i++)
      {
        Vector3D matCoord = phyxels[neighbour.get(i)].matCoord;
        vertex(matCoord.x, matCoord.y, matCoord.z);
      }
      
      
    }
    
    public void drawVertNeighbours(HE_Vertex vert)
    {
      Vector3D pos = new Vector3D (vert.xf(), vert.yf(), vert.zf());
      Phyxel chosenPhyxel = FindPhyxelAtLocation(pos);
      ArrayList< Integer > neighbour = GetNeighbours(chosenPhyxel,10);
      
      int iterNum = neighbour.size();
      strokeWeight(5);
      stroke(255, 0, 0);
      for(int i=0; i < iterNum; i++)
      {
        Vector3D matCoord = phyxels[neighbour.get(i)].matCoord;
        vertex(matCoord.x, matCoord.y, matCoord.z);
      }
    }
    
  
  }