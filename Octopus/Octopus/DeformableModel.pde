// A world consist of a spatial hash map //<>//
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
  
  int forceTime = 0;
  
  class DeformableModel {
    
    // Physics Parameters
    float Poisson_ratio;  // Poission Ratio
    float Young_modulus;  // Young's Modulus
    float massDensity;
    float massScaling;
    float damping = 500.0f;
    Vector3D C;           // Calculus combining the two above
    
    float kv;  // spring constant   ( Energy = 1/2 kx^2, x = delta_distance)
    
    
    // Spatial Hashing Map
    HashMap<Integer, ArrayList<Integer>> hashMap;
    public Phyxel[] phyxels;

    // World variables
    Vector3D worldSize;
    Vector3DI worldSizeInCell;
    float cellSize;
    float searchR;                  // search radius
    
    // Calculation Constants
    ArrayList<ArrayList<Float>> wlist;  // w_ij = W(|x_j - x_i|, h_i), W(r, h) is the polynomial kernel
    ArrayList<ArrayList<Vector3D>> xlist;  // x_ij = x_j - x_i
    WB_M33[] AInvlist;  // A_i = sum_j(x_ij * x_ij^T * w_ij) //<>//
    WB_Vector[] di;
    final int NEIGHBOURCOUNT = 10;
    float startTime = 0;
    final float DELTATIME = 0.002;
    float accumulator = 0;
    
    
    // Using HE_Mesh because it has no repeat vertex
    // Force the origin to be (0, 0, 0) and put the model in the positive areas
    public DeformableModel(Vector3D _worldSize, float _searchRadius, HE_Mesh _model, float _Young_modulus, float _Poisson_ratio)
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
      Young_modulus = _Young_modulus;
      Poisson_ratio = _Poisson_ratio;
      _loadModel(_model);
    }
    
    // Load the model into phyxels. 
    void _loadModel(HE_Mesh model)
    {
      WB_Coord [] vertices = model.getVerticesAsArray();
      //List<WB_Coord> vertexes = model.getPoints();
      int vertCount = vertices.length;
      print("vertCount", vertCount);
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
      
      // Calculate all the phyxel constants
      _initConstants();
    }
    
    void Update()
    {
      accumulator += (millis() - startTime) / 1000; 
      startTime = millis();

      if ( accumulator >= DELTATIME )
      {
        UpdatePhysics(DELTATIME);
        accumulator -= DELTATIME;
      }
      //startTime = millis();
      //print("deltaTime", deltaTime);
      //if (forceTime >= 5) return;
      
    }
    
    void UpdatePhysics(float deltaTime)
    {
      UpdateForces();
      Integrate(deltaTime);
    }
    
    void UpdateForces()
    {
      _initF();
      _updateF();
      printVec3D("f", phyxels[0].f);
      printVec3D("V", phyxels[0].v);
    }
    
    void Integrate(float deltaTime)
    {
      for( int i=0; i< phyxels.length; i++)
      {
        //printVec3D("f", phyxels[i].f);
        if (phyxels[i].f.magnitude() < 0.00001){
          phyxels[i].f = new Vector3D(0,0,0);
          phyxels[i].acc = new Vector3D(0,0,0);
          continue;
        }
        phyxels[i].acc = phyxels[i].f.mult(1/(phyxels[i].mass));
        Vector3D _a = phyxels[i].acc.mult(0.5*deltaTime*deltaTime);
        Vector3D _b = phyxels[i].v.mult(deltaTime);
        Vector3D _ab = (_b).add(_a);
        //phyxels[i].u = phyxels[i].u.add(_ab);
        phyxels[i].u = phyxels[i].u.add(_ab);
      }
      printVec3D("acc", phyxels[0].acc);
      
      UpdateForces();
      for(int i = 0; i < phyxels.length; i++)
      {
        if (phyxels[i].f.magnitude() < 0.00001){
          phyxels[i].f = new Vector3D(0,0,0);
          phyxels[i].acc = new Vector3D(0,0,0);
          continue;
        }
        Vector3D result = phyxels[i].f.mult(1/phyxels[i].mass);
        result = result.add(phyxels[i].acc);
        result = result.mult(deltaTime);
        result = result.mult(0.5);
        phyxels[i].v = result.add(phyxels[i].v);
      }
    }
    
    void _initF()
    {
      for (int i = 0; i < phyxels.length; i++)
      {
        phyxels[i].f = new Vector3D(0,0,0);
        phyxels[i].f = phyxels[i].f.subtract( phyxels[i].v.mult(damping));
      }
      //printVec3D("V", phyxels[0].v);
      printVec3D("finit", phyxels[0].f);
      
      // Add a continuous test force on one point
      //if(forceTime < 16){
        phyxels[0].f = new Vector3D(0, -100f, 0);
        forceTime ++;
      //}
    }
    
    /**************************************************************
    *                                                             *
    *                                                             *
    *                        HashMap                              *
    *                                                             *
    *                                                             *
    ***************************************************************/
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
    
    // Add a phyxel to all the cells it belongs to
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
    
    // Get the list of elements in the cell
    public ArrayList<Integer> GetListAt(int cell)
    {
      return hashMap.get(cell);
    }
    
    //public void ClearAt(int cell)
    //{
    //  hashMap.put (cell, new ArrayList<Integer>());
    //}
    
    // Clear all the elements in the hashMap
    public void ClearAll()
    {
      hashMap.clear();
      int total = worldSizeInCell.total();
      for(int i=0; i<total; i++)
      {
        hashMap.put (i, new ArrayList<Integer>());
      }
    }
    
    // Use a location to find the phyxel
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
    
    /**************************************************************
    *                                                             *
    *                                                             *
    *                          MATH                               *
    *                                                             *
    *                                                             *
    ***************************************************************/
    
    void _initConstants()
    {
      _initPhysicsParams();     
      _initNeighboursAndH();
      _initXW();
      _computeScalingFactor();
      _initMass();
      _initDensityVolume();
      _initADi();
    }
    
    // Calculate physics parameters
    void _initPhysicsParams()
    {
      //P = 0;        // The Poisson's ratio of a stable, isotropic, linear elastic material will be 
                       // greater than −1.0 or less than 0.5 because of the requirement for Young's modulus
      //E = pow(10, 4);
      massDensity = 1000.0f;
      kv = 100.0f;

      
      float d15 = Young_modulus / (1.0f + Poisson_ratio) / (1.0f - 2 * Poisson_ratio);
      float d16 = (1.0f - Poisson_ratio) * d15;
      float d17 = Poisson_ratio * d15;
      float d18 = Young_modulus / 2 / (1.0f + Poisson_ratio);
      
      C = new Vector3D(d16, d17, d18);
      
    }
    
    // Calculate h of each phyxel
    void _initNeighboursAndH()
    {
      for (int i = 0; i < phyxels.length; i++)
      {
        ArrayList<Integer> neighbours = GetNeighbours(phyxels[i], NEIGHBOURCOUNT);
        phyxels[i].setNeighbors(neighbours);
        
        float average_r = 0;
        for (int ii = 0; ii < neighbours.size(); ii++)
        {
          int j = neighbours.get(ii);
          average_r += (phyxels[i].matCoord).distance(phyxels[j].matCoord);
        }
        average_r /= neighbours.size();
        float h = average_r * 3;
        phyxels[i].h = h;
        phyxels[i].r = average_r;
      }
    }
    
    // Fill w_ij list and x_ij list and density and volume
    void _initXW() 
    { 
      wlist = new ArrayList<ArrayList<Float>>();
      xlist = new ArrayList<ArrayList<Vector3D>>();
      
      for (int i = 0; i<phyxels.length; i++) 
      {
        
        ArrayList<Float> _wlist = new ArrayList<Float>();
        ArrayList<Vector3D> _xlist = new ArrayList<Vector3D>();
        
        float h = phyxels[i].h;
        float h2 = h * h;
        float factor = 315.0f / (64 * (float)Math.PI * pow(h,9));
        
        ArrayList<Integer> neighbours = phyxels[i].getNeighbours();
        for(int ii = 0; ii < phyxels.length; ii++)
        {
          _wlist.add(0f);
          _xlist.add(new Vector3D(0,0,0));
        }
        
        int neighbourNum = neighbours.size();
        for (int ii = 0; ii < neighbourNum; ii++)
        {
          int j = neighbours.get(ii);
          float w = 0;
          Vector3D x_ij = (phyxels[j].matCoord).subtract(phyxels[i].matCoord); //  x_ij = x_j - x_i
          _xlist.set(j, x_ij);
          
          float r2 = x_ij.magnitude2();  // r = |x_j - x_i|
          if (r2 >= h2) w = 0;
          else w = factor *  pow((h2 - r2), 3);
          _wlist.set(j, w);
          
        }
        

        
        wlist.add(_wlist);
        xlist.add(_xlist);
      }
      
      //printArrayList("wlist", wlist.get(0));
    }
    
    // :)
    void _computeScalingFactor()
    {
      massScaling = 0.0f;
      int iterNum = phyxels.length;
      for(int i=0; i < iterNum; i ++)
      {
        float sum = 0;
        ArrayList<Integer>neighbours = phyxels[i].getNeighbours();
        int neighbourNum = neighbours.size();
        for(int ii=0; ii < neighbourNum; ii++)
        {
          int j = neighbours.get(ii);
          sum += pow(xlist.get(i).get(j).magnitude(), 3) * wlist.get(i).get(j);
        }
        massScaling += 1.0f / sum;
      }
      massScaling /= phyxels.length;
    }
    
    // :)
    void _initMass()
    {
      for(int i=0; i<phyxels.length; i++)
      {
        phyxels[i].mass = massScaling * pow(phyxels[i].r, 3) * massDensity;
      }
    }
    
    // :)
    void _initDensityVolume()
    {
      for(int i=0; i<phyxels.length; i++){
        float rho = 0;
        ArrayList<Integer> neighbours = phyxels[i].getNeighbours();
        int neighbourNum = neighbours.size();
        for(int ii=0; ii<neighbourNum; ii++)
        {
          int j = neighbours.get(ii);
          rho += phyxels[j].mass * wlist.get(i).get(j);
        }
        phyxels[i].density = rho;
        //print("density", i, rho, "\n");
        phyxels[i].volume = phyxels[i].mass / phyxels[i].density;
      }
      
    }
    
    // Calculate matrix A
    // A_i = sum_j(x_ij * x_ij^T * w_ij)
    void _initADi()
    {
      AInvlist = new WB_M33 [phyxels.length];
      di = new WB_Vector[phyxels.length];
      
      for(int i=0; i < phyxels.length; i++)
      {
        ArrayList<Integer> neighbours = phyxels[i].getNeighbours();
        int neighbourNum = neighbours.size();
        
        WB_M33 A = new WB_M33(0,0,0,0,0,0,0,0,0);
        WB_M33 a;

        for(int ii=0; ii < neighbourNum; ii++)
        {      
          int j = neighbours.get(ii);
          a = new WB_M33(xlist.get(i).get(j).x * xlist.get(i).get(j).x, 
                         xlist.get(i).get(j).x * xlist.get(i).get(j).y,
                         xlist.get(i).get(j).x * xlist.get(i).get(j).z,
                         xlist.get(i).get(j).y * xlist.get(i).get(j).x, 
                         xlist.get(i).get(j).y * xlist.get(i).get(j).y,
                         xlist.get(i).get(j).y * xlist.get(i).get(j).z,
                         xlist.get(i).get(j).z * xlist.get(i).get(j).x, 
                         xlist.get(i).get(j).z * xlist.get(i).get(j).y,
                         xlist.get(i).get(j).z * xlist.get(i).get(j).z);
          //printM33("a", a);
          a.mul(wlist.get(i).get(j));
          //printM33("a after multiply", a);
          A.add(a);
        }
        if(A.det() != 0){
          AInvlist[i] = A.inverse();
        }
        else{
          AInvlist[i] = new WB_M33(0,0,0,0,0,0,0,0,0);
          print("WARNING!!!!", i, "'s A has determinant 0!\n");
        }
        
        phyxels[i].dj = new ArrayList<WB_Vector>();
        di[i] = new WB_Vector(0, 0, 0);
        for(int ii=0; ii<neighbourNum; ii++)
        {
          
          int j = neighbours.get(ii);
          Vector3D xij_wij = (xlist.get(i).get(j)).mult(wlist.get(i).get(j));
          WB_Vector result = WB_M33.mulToVector(AInvlist[i], xij_wij.to_WB_Vector());
          phyxels[i].dj.add(result);
          di[i] = di[i].sub(result);
        }
      }
    }
    /*
    void _initM()
    {
      for (int i = 0; i < phyxels.length; i++)
      {
        phyxels[i].mass = 1;
      }
    }
    */
    
    
    
    //  Force of each phyxel has to be initialized to 0 at each frame
    //  (which has not been done yet!!!)
    void _updateF()
    {
      Vector3D _u, _v, _w;
      Vector3D du, dv, dw;
      WB_M33 J, JT;  // deriv U; deriv U transpose; Jacob
      WB_M33 epsilon = new WB_M33();
      WB_M33 sigma = new WB_M33();  // ε, σ
      WB_M33 Fe = new WB_M33();
      WB_M33 Fv = new WB_M33();
      WB_M33 Fe_add_Fv = new WB_M33();
      WB_Vector Ju, Jv, Jw;
      WB_Vector JvJw, JwJu, JuJv;
      //WB_Vector di, dj;
      WB_Vector _fi = new WB_Vector();
      WB_Vector _fj = new WB_Vector();
      for (int i=0; i < phyxels.length; i++)
      {
        _u = new Vector3D(0,0,0);
        _v = new Vector3D(0,0,0);
        _w = new Vector3D(0,0,0);
        ArrayList<Integer> neighbours = phyxels[i].getNeighbours();

        
        //_di = new Vector3D(0,0,0);
        // _di seems to be a constant too
        for (int ii=0; ii < neighbours.size(); ii++)
        {
          int j = neighbours.get(ii);
          _u = _u.add(xlist.get(i).get(j).mult((phyxels[j].u.x - phyxels[i].u.x) * wlist.get(i).get(j)));
          _v = _v.add(xlist.get(i).get(j).mult((phyxels[j].u.y - phyxels[i].u.y) * wlist.get(i).get(j)));
          _w = _w.add(xlist.get(i).get(j).mult((phyxels[j].u.z - phyxels[i].u.z) * wlist.get(i).get(j)));
          //_di = _di.add((xlist.get(i).get(j)).mult(-wlist.get(i).get(j)));
          //printVec3D("xlist_ij", xlist.get(i).get(j));
          //print("wlist_ij: " + wlist.get(i).get(j) + "\n");
          //print("i: " + i + ", j: " + ii + "\n"); 

        }
        
        /*
        if (i == 100){
          print("100 " +  phyxels[100].u.y + "\n");
          print("80 " +  phyxels[80].u.y + "\n");
          print("wlist_ij: " + wlist.get(100).get(80) + "\n");
          printVec3D ("mult", xlist.get(i).get(80).mult((phyxels[80].u.y - phyxels[100].u.y)));
          printVec3D ("calculation", xlist.get(i).get(80).mult((phyxels[80].u.y - phyxels[i].u.y) * wlist.get(i).get(80)));
          printVec3D ("_v", _v);
        }
        */

          

        
        //printVec3D("_u", _u);
        //printM33("Alist[i]", Alist.get(i));
        //WB_M33 A_inverse = Alist.get(i).inverse();
        WB_M33 A_inverse = AInvlist[i];
        
        /*
        if(i == 1){
          printVec3D("_u", _u);
          printVec3D("_v", _v);

          printVec3D("_w", _w);
          printM33("AInv", A_inverse);

        }
        */
        
        du = new Vector3D(A_inverse.m11 * _u.x + A_inverse.m12 * _u.y + A_inverse.m13 * _u.z,        
                          A_inverse.m21 * _u.x + A_inverse.m22 * _u.y + A_inverse.m23 * _u.z,
                          A_inverse.m31 * _u.x + A_inverse.m32 * _u.y + A_inverse.m33 * _u.z);
        dv = new Vector3D(A_inverse.m11 * _v.x + A_inverse.m12 * _v.y + A_inverse.m13 * _v.z,        
                          A_inverse.m21 * _v.x + A_inverse.m22 * _v.y + A_inverse.m23 * _v.z,
                          A_inverse.m31 * _v.x + A_inverse.m32 * _v.y + A_inverse.m33 * _v.z);
        dw = new Vector3D(A_inverse.m11 * _w.x + A_inverse.m12 * _w.y + A_inverse.m13 * _w.z,        
                          A_inverse.m21 * _w.x + A_inverse.m22 * _w.y + A_inverse.m23 * _w.z,
                          A_inverse.m31 * _w.x + A_inverse.m32 * _w.y + A_inverse.m33 * _w.z);
        /*
        du = new Vector3D(Alist.get(i).m11 * _u.x + Alist.get(i).m12 * _u.y + Alist.get(i).m13 * _u.z,        
                          Alist.get(i).m21 * _u.x + Alist.get(i).m22 * _u.y + Alist.get(i).m23 * _u.z,
                          Alist.get(i).m31 * _u.x + Alist.get(i).m32 * _u.y + Alist.get(i).m33 * _u.z);
        dv = new Vector3D(Alist.get(i).m11 * _v.x + Alist.get(i).m12 * _v.y + Alist.get(i).m13 * _v.z,        
                          Alist.get(i).m21 * _v.x + Alist.get(i).m22 * _v.y + Alist.get(i).m23 * _v.z,
                          Alist.get(i).m31 * _v.x + Alist.get(i).m32 * _v.y + Alist.get(i).m33 * _v.z);
        dw = new Vector3D(Alist.get(i).m11 * _w.x + Alist.get(i).m12 * _w.y + Alist.get(i).m13 * _w.z,        
                          Alist.get(i).m21 * _w.x + Alist.get(i).m22 * _w.y + Alist.get(i).m23 * _w.z,
                          Alist.get(i).m31 * _w.x + Alist.get(i).m32 * _w.y + Alist.get(i).m33 * _w.z);
        */
        
        J = new WB_M33(du.x + 1, du.y, du.z, dv.x, dv.y + 1, dv.z, dw.x, dw.y, dw.z + 1);
        /*
        if(i==1){
          printM33("J", J);
        }
        */
        
        JT = new WB_M33(0,0,0,0,0,0,0,0,0);
        J.transposeInto(JT);
        /*
        if(i == 1){
          printM33("JT", JT);
        }
        */
        epsilon = new WB_M33(0,0,0,0,0,0,0,0,0);
        epsilon = JT;
        epsilon.mul(J);
        //epsilon = WB_M33.mul(JT, J);
        epsilon.sub(new WB_M33(1,0,0,0,1,0,0,0,1));

        // Deriv u
        //printM33("U", U);
        // u transpose
        //UT = U;
        //UT.transpose();
        //printM33("UT", UT);
        // Jacobi
        
        // U * UT
        //UUT = WB_M33.mul(U, UT);
        
        
        // U = U + UT
        //U.add(UT);
        //epsilon = U + UT + U * UT
        //U.addInto(UUT, epsilon);
        //printM33("U", U);
        //printM33("UT", UT);
        //printM33("UUT", UUT);
        //printM33("epsilon", epsilon);
        // sigma = C * epsilon
        //sigma = epsilon;


        // Apply Calculus C
        //sigma.mul(C);
        
        sigma = new WB_M33(0,0,0,0,0,0,0,0,0);
        sigma = new WB_M33(
          C.x * epsilon.m11 + C.y * epsilon.m22 + C.y * epsilon.m33, C.z * epsilon.m12,                                     0,                                      
          0,                                                   C.y * epsilon.m11 + C.x * epsilon.m22 + C.y * epsilon.m33,  C.z * epsilon.m23, 
          C.z * epsilon.m31,                                   0,                                                          C.y * epsilon.m11 + C.y * epsilon.m22 + C.x * epsilon.m33
        );
        sigma.m13 = sigma.m31;
        sigma.m21 = sigma.m12;
        sigma.m32 = sigma.m23;
        
        /*
        sigma = epsilon;
        sigma = new WB_M33(
          sigma.m11 - Poisson_ratio * (sigma.m22 + sigma.m33), sigma.m12 * (1+Poisson_ratio),                       sigma.m13 * (1+Poisson_ratio),                                      
          sigma.m21 * (1+Poisson_ratio),                       sigma.m22 - Poisson_ratio * (sigma.m11 + sigma.m33), sigma.m23 * (1+Poisson_ratio), 
          sigma.m31 * (1+Poisson_ratio),                       sigma.m32 * (1+Poisson_ratio),                       sigma.m33 - Poisson_ratio * (sigma.m11 + sigma.m22)
        );
        sigma.mul(1/Young_modulus);
        */
        
        /*
        if(i==1){
          printM33("epsilon", epsilon);
          printM33("sigma", sigma);
        }
        */
        

        //printM33("sigma", sigma);
        Fe = new WB_M33(0,0,0,0,0,0,0,0,0);
        Fv = new WB_M33(0,0,0,0,0,0,0,0,0);

        WB_M33.mulInto(J, sigma, Fe);
        Fe.mul(-2 * phyxels[i].volume);
        /*
        if(i == 1){
          print("volume", phyxels[i].volume);
        }
        */
        //printM33("sigma", sigma);
        Ju = new WB_Vector(0,0,0);
        Jv = new WB_Vector(0,0,0);
        Jw = new WB_Vector(0,0,0);
        
        Ju = J.row(0);
        Jv = J.row(1);
        Jw = J.row(2);
        
        JvJw = new WB_Vector(0,0,0);
        JwJu = new WB_Vector(0,0,0);
        JuJv = new WB_Vector(0,0,0);
        
        
        JvJw = Jv.cross(Jw);
        JwJu = Jw.cross(Ju);
        JuJv = Ju.cross(Jv);
        
        /*
        if(i==1){
          printM33("J",J);
          printVecF("JvJw", JvJw);
          printVecF("JwJu", JwJu);
          printVecF("JuJv", JuJv);
          print("J.det", J.det());
        }
        */
        //di = WB_M33.mulToVector(Alist.get(i).inverse(), _di.to_WB_Vector());  //  di = A^-1*(-sum_j(x_ij*w_ij))
        //printVecF("di", di);
        
        

        //printM33("J", J);
        //printM33("Fe", Fe);
        //printM33("Fe", Fe);
        Fv = new WB_M33(JvJw.coords()[0], JvJw.coords()[1], JvJw.coords()[2],
                        JwJu.coords()[0], JwJu.coords()[1], JwJu.coords()[2],
                        JuJv.coords()[0], JuJv.coords()[1], JuJv.coords()[2]);
        //Fv.transpose();
        
        Fv.mul(-phyxels[i].volume * kv * (J.det() - 1));
        
       
        
        Fe.addInto(Fv, Fe_add_Fv);
        //printM33 ("Fe", Fe);
        //printM33 ("Fe+Fv", Fe_add_Fv);
        _fi = new WB_Vector(0,0,0);
        _fi = WB_M33.mulToVector(Fe_add_Fv, di[i]);
        //_fi = WB_M33.mulToVector(Fv,di);
        
        /*
        if(i == 1){
          printM33("Fe", Fe);
          printM33("Fv", Fv);
          printVec3D("u", phyxels[i].u);
          printVecF("di", di[i]);
        }
        */
        
        //printVecF ("_fi", _fi);
        //printVec3D ("f", phyxels[i].f);
        phyxels[i].f = phyxels[i].f.add(new Vector3D(_fi.coords()));
        /*
        if(i==1){
          printVec3D ("fi", phyxels[i].f);
        }
        */
        for (int ii=0; ii < neighbours.size(); ii++)
        {
          int j = neighbours.get(ii);
          //dj = WB_M33.mulToVector(Alist.get(i).inverse(), (xlist.get(i).get(j)).mult(wlist.get(i).get(j)).to_WB_Vector());
          //_fj = WB_M33.mulToVector(Fe_add_Fv, dj);
          //_fj = WB_M33.mulToVector(Fe_add_Fv, dj);
          _fj = WB_M33.mulToVector(Fe_add_Fv, phyxels[i].dj.get(ii));
          phyxels[j].f = phyxels[j].f.add(new Vector3D(_fj.coords()));
        }
        //printVec3D ("f", phyxels[i].f);
        
      }
      
    }
    

    /**************************************************************
    *                                                             *
    *                                                             *
    *                        VISUALIZATION                        *
    *                                                             *
    *                                                             *
    ***************************************************************/
    
    // Display the nearest neighbours on the static model
    public void draw(int idx)
    {
      int iterNum = phyxels.length;
      stroke(0, 0, 0);
      for(int i=0; i < iterNum; i++)
      {
        Vector3D matCoord = phyxels[i].matCoord;
        Vector3D u = phyxels[i].u;
        vertex(matCoord.x + u.x , matCoord.y + u.y , matCoord.z + u.z );
      }
      
      strokeWeight(10);
      stroke(0,255,0);
      vertex(phyxels[idx].matCoord.x, phyxels[idx].matCoord.y, phyxels[idx].matCoord.z);
      ArrayList< Integer > neighbour = GetNeighbours(phyxels[idx], NEIGHBOURCOUNT);
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
      ArrayList<Integer> neighbour = GetNeighbours(chosenPhyxel,10);
      
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