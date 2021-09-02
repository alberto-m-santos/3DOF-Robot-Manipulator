const
 NumJoints = 3;
 // Duration of the speed simulation
 Duration_speed_simulation = 10;
 // Initial angles of the joints
 initial_angle_joint0 = 45;
 initial_angle_joint1 = 30;
 initial_angle_joint2 = -40;
 // Tolerance of detecting singularity
 singularity_tolerance = 4e-2;

// Global Variables
var
  irobot, iB4: integer;
  d1, d2, d3: double;
  
  R01, R12, R23: Matrix;
  R03: Matrix;

  ReqThetas: matrix;
  JointsVelocities: matrix;
  Tol, tis: double;
  elapsed_time_speed_simulation: double;
  EndEffectorVelocity: matrix;
  singularity: boolean;

// Computes the cross product between two column matrixes
function crossProduct(Matrix1, Matrix2: matrix): matrix;
var
a,b,c,d,e,f: Double;
crossProductMatrix: matrix;
begin
a := Mgetv(Matrix1,0,0);
b := Mgetv(Matrix1,1,0);
c := Mgetv(Matrix1,2,0);
d := Mgetv(Matrix2,0,0);
e := Mgetv(Matrix2,1,0);
f := Mgetv(Matrix2,2,0);

crossProductMatrix := Mzeros(3,1);
Msetv(crossProductMatrix, 0, 0, b*f - e*c);
Msetv(crossProductMatrix, 1, 0, c*d - a*f);
Msetv(crossProductMatrix, 2, 0, a*e - b*d);

result := crossProductMatrix;
end;

// Computes the determinant of a matrix 3x3
function Mdet3x3(mat: matrix): double;
var
  a11,a12,a13,a21,a22,a23,a31,a32,a33: double;
begin
  a11 := Mgetv(mat, 0,0);
  a12 := Mgetv(mat, 0,1);
  a13 := Mgetv(mat, 0,2);

  a21 := Mgetv(mat, 1,0);
  a22 := Mgetv(mat, 1,1);
  a23 := Mgetv(mat, 1,2);

  a31 := Mgetv(mat, 2,0);
  a32 := Mgetv(mat, 2,1);
  a33 := Mgetv(mat, 2,2);

  result := (a11*a22*a33 + a12*a23*a31 + a13*a21*a32) - (a31*a22*a13 + a32*a23*a11 + a33*a21*a12);
end;


function DHMat(a, alpha, d, theta: double): Matrix;
var ct, st, ca, sa: double;
    R: Matrix;
begin
  ct := cos(theta);
  st := sin(theta);
  ca := cos(alpha);
  sa := sin(alpha);
  
  R := Meye(4);
  MSetV(R, 0, 0, ct); MSetV(R, 0, 1,-st * ca);  MSetV(R, 0, 2, st * sa);  MSetV(R, 0, 3, a * ct);
  MSetV(R, 1, 0, st); MSetV(R, 1, 1, ct * ca);  MSetV(R, 1, 2,-ct * sa);  MSetV(R, 1, 3, a * st);
  MSetV(R, 2, 0,  0); MSetV(R, 2, 1, sa     );  MSetV(R, 2, 2, ca     );  MSetV(R, 2, 3, d     );

  result := R;
end;


function RotZMat(theta: double): Matrix;
var ct, st: double;
    R: Matrix;
begin
  ct := cos(theta);
  st := sin(theta);

  R := Meye(3);
  MSetV(R, 0, 0, ct); MSetV(R, 0, 1,-st);  MSetV(R, 0, 2, 0);
  MSetV(R, 1, 0, st); MSetV(R, 1, 1, ct);  MSetV(R, 1, 2, 0);
  MSetV(R, 2, 0,  0); MSetV(R, 2, 1, 0 );  MSetV(R, 2, 2, 1);

  result := R;
end;


function RotXMat(theta: double): Matrix;
var ct, st: double;
    R: Matrix;
begin
  ct := cos(theta);
  st := sin(theta);

  R := Meye(3);
  MSetV(R, 0, 0, 1 ); MSetV(R, 0, 1, 0 );  MSetV(R, 0, 2, 0  );
  MSetV(R, 1, 0, 0 ); MSetV(R, 1, 1, ct);  MSetV(R, 1, 2, -st);
  MSetV(R, 2, 0,  0); MSetV(R, 2, 1, st);  MSetV(R, 2, 2, ct );

  result := R;
end;

// Place here the Inverse Kinematics calculations
// Return a column matrix with Theta1, Theta2 and Theta3
function IK3(XYZ: matrix): matrix;
var
  beta, phi, s, r, c, d, x, y, z, num, den, Theta1, Theta2, Theta3: double;
  Thetas: Matrix;
begin
  // desired positions
  x := Mgetv(XYZ, 0, 0);
  y := Mgetv(XYZ, 1, 0);
  z := Mgetv(XYZ, 2, 0);

  // auxiliar terms
  s := z - d1;
  r := sqrt(sqr(x) + sqr(y));
  c := sqrt(sqr(s) + sqr(r));
  num := sqr(c) - Sqr(d2) - Sqr(d3);
  den := 2 * d2 * d3;
  phi := Atan2(s,r);
  beta := arccos((sqr(c) + Sqr(d2) - Sqr(d3)) / (2 * d2 * c));

  // angles computation
  Theta1 := Atan2(y,x);
  Theta2 := -(phi - beta);
  Theta3 := -arccos(num/den);
  Thetas := Mzeros(3, 1);
  Msetv(Thetas, 0, 0, Theta1);
  Msetv(Thetas, 1, 0, Theta2);
  Msetv(Thetas, 2, 0, Theta3);

  //return the angles vector
  result := Thetas;
end;


// Place here the Direct Kinematics calculations
// Return a column matrix with X, X, and Z
function DK3(Thetas: matrix): matrix;
var
  A01, A12, A23, A03 :Matrix; // DH matrix variables
  XYZ_column_vector :Matrix;  // position of the end effector vector
begin
  // Compute the homogeneous DH transformation matrixes
  // using DHMat function:
  // parameters: a, alpha, d, theta
  // return: DH transformation matrix
  A01 := DHMat(0, -pi*0.5, d1 , Mgetv(Thetas,0,0));
  A12 := DHMat(d2, 0, 0, Mgetv(Thetas,1,0));
  A23 := DHMat(d3, 0, 0, Mgetv(Thetas,2,0));

  // Compute the DH transformation matrix between base joint and end-effector
  A03 := MMult(MMult(A01,A12),A23);

  // Get the position of the end-effector from the last column of A03 matrix
  XYZ_column_vector := Mzeros(3,1);
  Msetv(XYZ_column_vector, 0, 0, Mgetv(A03, 0, 3));
  Msetv(XYZ_column_vector, 1, 0, Mgetv(A03, 1, 3));
  Msetv(XYZ_column_vector, 2, 0, Mgetv(A03, 2, 3));

  // return the position vector
  result := XYZ_column_vector;
end;

procedure printMatrix(m: matrix; row,col :integer);
begin
  MatrixToRange(row,col,m);
end;

procedure printValue(value: double; row,col :integer);
begin
  SetRCValue(row,col,format('%.3g',[value]));
end;

procedure printText(text: string; row,col :integer);
begin
  SetRCValue(row,col,text);
end;
function computeJacobian: matrix;
var
  Theta1, Theta2, Theta3: double;
  A01, A12, A23, A02 ,A03: matrix;
  R00,R01,R02 :matrix;
  d00, d01, d02, d03: matrix;
  ZColumnVector: matrix;
  Jv1, Jv2, Jv3: matrix;
  A, B: matrix;
  Jacobian :matrix;
  invertedJacobian :matrix;
begin
  //Z unit vector
  ZColumnVector := Mzeros(3,1);
  Msetv(ZColumnVector, 2, 0, 1);

  //Current angles of the joints
  Theta1 := GetAxisPos(irobot, 0);
  Theta2 := GetAxisPos(irobot, 1);
  Theta3 := GetAxisPos(irobot, 2);

  //Compute the Homogeneous Transformation Matrixes using DH method
  A01 := DHMat(0, -pi*0.5, d1 , Theta1);
  A12 := DHMat(d2, 0, 0, Theta2);
  A23 := DHMat(d3, 0, 0, Theta3);

  //Jacobian columns
  //Jv1
  R00 := Meye(3);
  A02 := MMult(A01,A12);
  A03 := MMult(A02,A23);
  d03 := MCrop(A03, 0, 3, 2, 3);
  d00 := Mzeros(3,1);
  A := MMult(R00, ZColumnVector);
  B := Msub(d03, d00);
  Jv1 := crossProduct(A, B);

  //Jv2
  R01 := Mcrop(A01, 0, 0, 2, 2);
  d01 := MCrop(A01, 0, 3, 2, 3);
  A := Mmult(R01, ZColumnVector);
  B := Msub(d03, d01);
  Jv2 := crossProduct(A, B);

  //Jv3
  R02 := MCrop(A02, 0, 0, 2, 2);
  d02 := MCrop(A02, 0, 3, 2, 3);
  A03 := Mmult(A02, A23);
  d03 := MCrop(A03, 0, 3, 2, 3);
  A := Mmult(R02, ZColumnVector);
  B := Msub(d03,d02);
  Jv3 := crossProduct(A, B);

  Jacobian := Mzeros(3,3);
  Msetv(Jacobian, 0, 0, Mgetv(Jv1, 0, 0));
  Msetv(Jacobian, 1, 0, Mgetv(Jv1, 1, 0));
  Msetv(Jacobian, 2, 0, Mgetv(Jv1, 2, 0));
  Msetv(Jacobian, 0, 1, Mgetv(Jv2, 0, 0));
  Msetv(Jacobian, 1, 1, Mgetv(Jv2, 1, 0));
  Msetv(Jacobian, 2, 1, Mgetv(Jv2, 2, 0));
  Msetv(Jacobian, 0, 2, Mgetv(Jv3, 0, 0));
  Msetv(Jacobian, 1, 2, Mgetv(Jv3, 1, 0));
  Msetv(Jacobian, 2, 2, Mgetv(Jv3, 2, 0));

  result := Jacobian;
end;

// Check singularity by looking at the jacobian
function checkSingularity(jacobian :matrix): boolean;
begin
if abs(Mdet3x3(Jacobian)) < singularity_tolerance then begin
    singularity := true;
    result := true;
  end
  else begin
    singularity := false;
    result := false;
  end;
end;

// Using the inverted jacobian computes and returns the joint velocities
// Takes as input the desired velocity of the end effector
function ComputeJointVelocities(endEffectorVelocity: matrix): matrix;
var
  Jacobian :matrix;
  invertedJacobian :matrix;
begin
   Jacobian := computeJacobian;
  //Check singularity
  if checkSingularity(Jacobian) then begin
    result := Mzeros(3,1);
  end
  else begin
    invertedJacobian := Minv(Jacobian);
    result := Mmult(invertedJacobian, endEffectorVelocity);
  end;
end;



procedure SetThetas(Thetas: matrix);
var i: integer;
begin
  for i := 0 to NumJoints - 1 do begin
    SetMotorControllerMode(iRobot, i,'pidposition');
    //                                ki    kd   kp  kf
    SetMotorControllerPars(iRobot, i, 0.01, 20, 100, 0);
    SetAxisPosRef(iRobot, i, Mgetv(Thetas, i, 0));
  end;
end;

procedure SetVelocity(Velocity: matrix);
var i: integer;
begin
  for i := 0 to NumJoints - 1 do begin
    SetMotorControllerMode(iRobot, i,'pidspeed');
    //                                ki    kd   kp  kf
    SetMotorControllerPars(iRobot, i, 0.01, 0, 4, 0);
    SetAxisSpeedRef(iRobot, i, Mgetv(Velocity, i, 0));
  end;
end;

function JointError(ReqThethas: matrix): double;
var err: double;
    i: integer;
begin
  err := 0;
  for i := 0 to NumJoints - 1 do begin
    err := err + abs(GetAxisPos(iRobot, i) - Mgetv(ReqThethas, i, 0));
  end;
  result := err;
end;

//Reset the joints angles to a predefined position under constant values named:
// initial_angle_joint0
// initial_angle_joint1
// initial_angle_joint2
procedure reset;
var initialThetas: matrix;
begin
  initialThetas := Mzeros(3,1);
  Msetv(initialThetas, 0, 0,rad(initial_angle_joint0));
  Msetv(initialThetas, 1, 0,rad(initial_angle_joint1));
  Msetv(initialThetas, 2, 0,rad(initial_angle_joint2));
  SetThetas(initialThetas);
end;

procedure Control;
var i: integer;
    xw, yw, zw: double;
    t1, t2, t3: double;
    XYZ: matrix;
    ReqThetasDeg: Matrix;
begin
  // Print end effector position to the sheet
  MatrixToRange(11, 2, GetSolidPosMat(iRobot, iB4));

  // Print end effector velocity to the sheet
  MatrixToRange(7,2, GetSolidLinearVelMat(iRobot, iB4));

  // Read joint positions and show
  for i := 0 to NumJoints -1 do begin
    //imprimir na 2a coluna e da 3a linha à 5a linha na sheet as
    // posiçoes angulares das joints
    SetRCValue(3 + i, 2, format('%.3g',[Deg(GetAxisPos(irobot, i))]));
  end;

  // SE SE PREMIR O BOTÃO Direct
  if RCButtonPressed(2, 4) then begin
    // Cria a matriz coluna ReqThetas como input para IK3
    ReqThetas := Mzeros(3, 1);
    Msetv(ReqThetas, 0, 0, rad(GetRCValue(3, 4)));
    Msetv(ReqThetas, 1, 0, rad(GetRCValue(4, 4)));
    Msetv(ReqThetas, 2, 0, rad(GetRCValue(5, 4)));

    // Move o robô para os angulos inseridos na sheet
    SetThetas(ReqThetas);

    // Faz a CINEMÁTICA DIRETA e retorna a matriz coluna XYZ
    XYZ := DK3(ReqThetas);

    // imprimir XYZ
    MatrixToRangeF(11, 3, XYZ, '%.3f');
  end;

  // SE SE PREMIR O BOTÃO indirect:
  if RCButtonPressed(10,4) then begin
      // Cria a matriz coluna XYZ como input para IK3
      XYZ := Mzeros(3,1);

      Msetv(XYZ, 0, 0, GetRCValue(11,4));
      Msetv(XYZ, 1, 0, GetRCValue(12,4));
      Msetv(XYZ, 2, 0, GetRCValue(13,4));

      // Faz a CINEMÁTICA INVERSA e retorna a matriz coluna com os Thetas
      ReqThetas := IK3(XYZ);

      // Move o robô para os angulos retornados por IK3
      SetThetas(ReqThetas);

      // imprime para a sheet o vetor ReqThetas em graus
      ReqThetasDeg := Mzeros(3,1);
      Msetv(ReqThetasDeg, 0, 0, Deg(Mgetv(ReqThetas, 0, 0)));
      Msetv(ReqThetasDeg, 1, 0, Deg(Mgetv(ReqThetas, 1, 0)));
      Msetv(ReqThetasDeg, 2, 0, Deg(Mgetv(ReqThetas, 2, 0)));
      MatrixToRangeF(3, 3, ReqThetasDeg, '%.1f');
  end;

  // Se se premir o botão SetVelocity:
  if RCButtonPressed(6,5) then begin
    // * desativa singularidade
    // * começa tempo de simulaçao
    singularity := false;
    JointsVelocities := Mzeros(3,1);
    elapsed_time_speed_simulation := 0;

    Msetv(EndEffectorVelocity, 0, 0, GetRCValue(7,5));
    Msetv(EndEffectorVelocity, 1, 0, GetRCValue(8,5));
    Msetv(EndEffectorVelocity, 2, 0, GetRCValue(9,5));

    // calcula as velocidades de articulaçoes
    // em funçao da velocidade do efetor solicitada pelo utilizador na sheet
    JointsVelocities := ComputeJointVelocities(EndEffectorVelocity);

    // aplica as velocidades no controlador
    setVelocity(JointsVelocities);

    // imprime as velocidades das articulaçoes calculadas para a sheet
    printMatrix(JointsVelocities,7,4);
  end;

  // O que fazer em caso de controlo de velocidade:
  if GetMotorControllerMode(irobot,0) = 'PIDSpeed' then begin
    // Verifica em tempo real se ocorre singularidade
    checkSingularity(computeJacobian);

    // Se o tempo de simulaçao acabou, pára o motor e volta à posiçao inicial
    // passando para controlo de posiçao outra vez
    if elapsed_time_speed_simulation >  Duration_speed_simulation then begin
      elapsed_time_speed_simulation := 0;
      reset;
    end
    else begin
    // Se acontecer singularidade entao pára motor e volta à posiçao inicial
      if singularity then begin
         elapsed_time_speed_simulation := 0;
         setVelocity(Mzeros(3,1));
         reset;
      end
      else begin
      // Se nao acontecer singularidade e se simuçao a decorrer, aplica as velocidades
      // nas articulaçoes. Incrementa tempo de simulaçao.
        setVelocity(JointsVelocities);
        printValue(elapsed_time_speed_simulation,3,9);
        elapsed_time_speed_simulation := elapsed_time_speed_simulation + 0.04;
      end;
    end;
  end;

  // Se acontecer singularidade imprime: sim.
  // Se não acontecer, imprime: não.
  if singularity then begin
    printText('YES',3,6);
  end
  else begin
    printText('NO',3,6);
  end;

  setRCValue(3,7,GetMotorControllerMode(irobot,0));
  printValue(tis,3,8);
  tis := tis + 0.04;
end;

procedure Initialize;
var i: integer;
begin
  irobot := 0;

  iB4 := GetSolidIndex(irobot, 'B4');

  SetRCValue(2, 1, 'Joint');
  SetRCValue(2, 2, 'Pos (deg)');
  for i := 0 to NumJoints -1 do begin
    SetRCValue(3 + i, 1, format('%d',[i]));
  end;

  // DK
  SetRCValue(10, 3, 'DK');

  // IK
  SetRCValue(2, 3, 'IK');

  // Velocity
  SetRCValue(6, 3, 'Vjoints');
  SetRCValue(7, 3, 'v1:');
  SetRCValue(8, 3, 'v2:');
  SetRCValue(9, 3, 'v3');

  SetRCValue(6, 1, 'Veffector');
  SetRCValue(7, 1, 'Vx:');
  SetRCValue(8, 1, 'Vy:');
  SetRCValue(9, 1, 'Vz:');


  SetRCValue(6, 5, '[SetVelocity]');

  printText('Singularity',2,6);
  printText('Working Mode',2,7);
  printText('Simulation time',2,8);
  printText('Jacobian simulation time',2,9);

  d1 := 0.55;
  d2 := 0.4;
  d3 := 0.325;

  Tol := 0.2;

  // Contador de tempo de simulaçao
  elapsed_time_speed_simulation := 0;

  // Vetor global das velocidades do efetuador
  EndEffectorVelocity := Mzeros(3,1);

  // Booleano sinalizador de singularidade
  singularity := false;
end;
