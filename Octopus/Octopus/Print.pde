Boolean printOn = true;

void printM33(String name, WB_M33 matrix)
{
  if(!printOn) return;
  print (name + ":\n");
  print(matrix.m11 + " " + matrix.m12 + " " + matrix.m13 + "\n");
  print(matrix.m21 + " " + matrix.m22 + " " + matrix.m23 + "\n");
  print(matrix.m31 + " " + matrix.m32 + " " + matrix.m33 + "\n\n");
}

void printVecF(String name, WB_Vector vec)
{
  if(!printOn) return;
  print (name + ":\n");
  print(vec.xf() + " " + vec.yf() + " " + vec.zf() + "\n\n");
 
}

void printVec3D(String name, Vector3D vec)
{
  if(!printOn) return;
  print (name + ":\n");
  print(vec.x + " " + vec.y + " " + vec.z + "\n\n");
}

void printVar(String name, String variable)
{
  if(!printOn) return;
  print(name + ":\n");
  print(variable +"\n");
}

void printVar(String name, int variable)
{
  if(!printOn) return;
  print(name + ":\n");
  print(variable + "\n");
}

void printVar(String name, float variable)
{
  if(!printOn) return;
  print(name + ":\n");
  print(variable + "\n");
}

void printVar(String name, double variable)
{
  if(!printOn) return;
  print(name + ":\n");
  print(variable + "\n");
}

void printArrayList(String name, ArrayList<Float> list)
{
  if(!printOn) return;
  print(name + ":\n");
  for(int i = 0; i < list.size() / 8; i++)
  {
    for(int j = 0; j < 8; j++){
      print(list.get(i*8+j));
      print(" ");
    }
    print("\n");
  }
}

void printArrayListI(String name, ArrayList<Integer> list)
{
  if(!printOn) return;
  print(name + ":\n");
  for(int i = 0; i < list.size() / 8; i++)
  {
    for(int j = 0; j < 8; j++){
      print(list.get(i*8+j));
      print(" ");
    }
    print("\n");
  }
}