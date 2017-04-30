void printM33(String name, WB_M33 matrix)
{
  print (name + ":\n");
  print(matrix.m11 + " " + matrix.m12 + " " + matrix.m13 + "\n");
  print(matrix.m21 + " " + matrix.m22 + " " + matrix.m23 + "\n");
  print(matrix.m31 + " " + matrix.m32 + " " + matrix.m33 + "\n\n");
}

void printVecF(String name, WB_Vector vec)
{
  print (name + ":\n");
  print(vec.xf() + " " + vec.yf() + " " + vec.zf() + "\n\n");
 
}

void printVec3D(String name, Vector3D vec)
{
  print (name + ":\n");
  print(vec.x + " " + vec.y + " " + vec.z + "\n\n");
 
}

void printVar(String name, String variable)
{
  print(name + ":\n");
  print(variable);
}