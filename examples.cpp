#include "Eigen/Dense"

#include <iostream>

int main()
{
    // NOTE: auto makes a reference, explicit type definition is a copy !!

    // base Matrix declaration

    // Matrices with size set at compile time (in stack)

    // Matrix<typename Scalar, int RowsAtCompileTime, int ColsAtCompileTime>
    Eigen::Matrix<double, 4, 4> double4x4matrixA;

    // same as typedef: Matrix4d
    Eigen::Matrix4d double4x4matrixB;

    // Vector3d is a column vector. It is actually a base Matrix
    // these 2 are the same
    Eigen::Vector3d vectorA;
    Eigen::Matrix<double, 3, 1> vectorB;

    // RowVector3d is a row vector. It is actually a base Matrix
    // these 2 are the same
    Eigen::RowVector3d rowVectorA;
    Eigen::Matrix<double, 1, 3> rowVectorB;

    // Matrices with size set during run time (in heap)
    // they use size variables set during run time
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> doubleDynMatrixA;

    // same as typedef: MatrixXd
    Eigen::MatrixX4d doubleDynMatrixB;

    // these 2 are also the same
    Eigen::Matrix<double, Eigen::Dynamic, 1> vectorDynA;
    Eigen::VectorXd vectorDynB;

    // default constructors
    // this one allocates space in the stack but does not initialize it
    Eigen::Matrix4d m4dA;
    // this one does not allocate dynamic space (heap) and does not initialize it
    Eigen::MatrixXd mXdA;

    // constructors with dynamic size. not initialized
    Eigen::MatrixXd mXdB(3,3);
    Eigen::VectorXd vXdA(4);

    // non dynamic construct with initialized values
    //Eigen::Matrix4d m4dB(1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1);
    Eigen::Vector4d v4dA {1,2,3,4};
    Eigen::Matrix4d m4dB {{1,1,1,1}, {2,2,2,2}, {3,3,3,3}, {4,4,4,4}};

    // if you group initialize a dynamic matrix it will set it's size correctly
    Eigen::MatrixXd mXdC {{1,2}, {3,4}, {5,6}};

    // ad dynamic vector needs to be group initialized
    Eigen::VectorXd vXdB {{1.5, 2.5, 3.5}};
    Eigen::RowVectorXd rvXdA {{2.5, 3.5, 4.5}};

    // assignement will copy size of right hand into left hand dynamic matrix
    Eigen::Matrix3f rhm;
    Eigen::MatrixXf lhm(2,2);
    std::cout << "\nassignement of matrix with different size matrix:\n";
    std::cout << "size = " << lhm.rows() << " " << lhm.cols() << std::endl;
    lhm = rhm;
    std::cout << "size = " << lhm.rows() << " " << lhm.cols() << std::endl;

    // predefined matrices
    std::cout << "\npredefined matrices:\n";
    std::cout << std::endl << Eigen::Matrix4f::Random() << std::endl;
    std::cout << std::endl << Eigen::Matrix4f::Ones() << std::endl;
    std::cout << std::endl << Eigen::Matrix4f::Zero() << std::endl;
    std::cout << std::endl << Eigen::Matrix4f::Constant(1.2) << std::endl;
    std::cout << std::endl << Eigen::RowVectorXi::LinSpaced(4,2,15) << std::endl;
    std::cout << std::endl << Eigen::RowVector4d::UnitX() << std::endl;
    std::cout << std::endl << Eigen::RowVector4d::UnitY() << std::endl;
    std::cout << std::endl << Eigen::RowVector4d::UnitZ() << std::endl;
    std::cout << std::endl << Eigen::RowVector4d::UnitW() << std::endl;
    std::cout << std::endl << Eigen::RowVector4d::Unit(1) << std::endl;
    
    // for non linear algebra we have arrays to perform coefficient wise operations
    Eigen::Array<double, 3, 4> double3x4arrayA;

    // these arrays are the same
    Eigen::Array<float,Eigen::Dynamic,1> arrayA;                Eigen::ArrayXf arrayB;
    Eigen::Array<float,3,1> arrayC;                             Eigen::Array3f arrayD;
    Eigen::Array<double,Eigen::Dynamic,Eigen::Dynamic> arrayE;  Eigen::ArrayXXd arrayF;
    Eigen::Array<double,3,3> arrayG;                            Eigen::Array33d arrayH;

    // example of coefficent wise multiplication instead of matrix mult
    Eigen::ArrayXXd aXdA {{1,2}, {3,4}, {5,6}};
    Eigen::ArrayXXd aXdB {{1,2}, {3,4}, {5,6}};
    auto aXdC = aXdA * aXdB;
    std::cout << "\nmultiplying arrays:\n";
    std::cout << aXdC << std::endl;

    // some other coefficient wise operations
    Eigen::Array<double, 2,3> arrayTest;
    arrayTest.setRandom(); // between -1 and 1
    std::cout << "\narray abs():\n";
    std::cout << arrayTest.abs() << std::endl;
    std::cout << "\narray.abs().sqrt():\n";
    std::cout << arrayTest.abs().sqrt() << std::endl;

    Eigen::ArrayXXd aXdD(3,2);
    aXdD.setConstant(3);
    std::cout << "\nget smalles values between arrays array.min(array):\n";
    std::cout << aXdA.min(aXdD) << std::endl;
    std::cout << "\nget biggest values between arrays array.min(array):\n";
    std::cout << aXdA.max(aXdD) << std::endl;

    // block operations. they work both in Matrix and Array classes

    Eigen::MatrixXf blockTest(4,4);
    blockTest <<
        1, 2, 3, 4,
        5, 6, 7, 8,
        9,10,11,12,
        13,14,15,16;
    std::cout << "\nblock in the middle:\n";
    Eigen::MatrixXf blockA = blockTest.block<2,2>(1,1); // if size is known during compile
    std::cout << blockA << std::endl;
    auto blockB = blockTest.block<2,2>(1,1); // if size is known during compile
    blockA(0,0) = 66; // this is a copy
    std::cout << blockTest << std::endl;
    blockB(0,0) = 77; // this is a reference
    std::cout << blockTest << std::endl;

    // getting a dynamic sized block is
    Eigen::MatrixXf blockC = blockTest.block(1,1,2,2); // 2 and 2 not known in advance

    // changing block inside an array
    Eigen::Array22f smallBlock;
    smallBlock << 1,2, 3,4;
    Eigen::Array44f largeBlock = Eigen::Array44f::Constant(9);
    std::cout << "\nblock change original:\n";
    std::cout << largeBlock << std::endl;
    std::cout << "\nblock change middle:\n";
    std::cout << smallBlock << std::endl;
    largeBlock.block<2,2>(1,1) = smallBlock;
    std::cout << "\nblock change result :\n";
    std::cout << largeBlock << std::endl;

    // this could be used in camera projection of vertices
    Eigen::Array<double, 4, 3> vertices;
    vertices.row(0).setConstant(8);
    vertices.row(1).setConstant(16);
    vertices.row(2).setConstant(32);
    vertices.row(3).setLinSpaced(3,1,3);
    vertices.row(3) = Eigen::Array3d::Constant(2).pow(Eigen::Array3d::LinSpaced(3,1,3));
    std::cout << "\nvertices\n";
    std::cout << vertices.transpose() << std::endl;
    for ( auto col : vertices.colwise() ) col /= col(3);
    std::cout << vertices.transpose() << std::endl;

    // some cool shit
    Eigen::Matrix3f m = Eigen::Matrix3f::Random();
    std::ptrdiff_t i, j;
    float minOfM = m.minCoeff(&i,&j);
    std::cout << "\nHere is the matrix m:\n" << m << std::endl;
    std::cout << "Its minimum coefficient (" << minOfM 
       << ") is at position (" << i << "," << j << ")\n\n";

    // Conversions between matrix and array world
    std::cout << "\nCoefficient wise multiplication between arrays\n";
    Eigen::Array<float, 2, 2> a1, a2, a3;
    a1 << 1, 2, 3, 4;
    a2 << 3, 4, 5, 6;
    std::cout << "\na1:\n" << a1 << std::endl;
    std::cout << "\na2:\n" << a2 << std::endl;
    a3 = a1 * a2;
    std::cout << "\na3:\n" << a3 << std::endl;

    std::cout << "\nMatrix product between matrices\n";
    Eigen::Matrix<float, 2, 2> m1, m2, m3;
    m1 << 1, 2, 3, 4;
    m2 << 3, 4, 5, 6;
    std::cout << "\nm1:\n" << m1 << std::endl;
    std::cout << "\nm2:\n" << m2 << std::endl;
    m3 = m1 * m2;
    std::cout << "\nm3:\n" << m3 << std::endl;

    // mixing array and matrix is forbidden
    //auto test = a1 + m2;

    // you can convert back and forth though
    a3 = a1 + m2.array();
    m3 = a1.matrix() + m2;

    // size constructor with Dynamic size, constructed in heap, not stack
    Eigen::Matrix<float, Eigen::Dynamic, 3>  dz1(2,3);
    dz1.setOnes(); // values have to be initialized
    std::cout << "\ndynamic size constructor A:\n" << dz1 << std::endl;
    Eigen::Matrix<float, 3, Eigen::Dynamic>  dz2(3,2);
    dz2.setOnes(); // values have to be initialized
    std::cout << "\ndynamic size constructor B:\n" << dz2 << std::endl;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>  dz3(4,5);
    dz3.setOnes(); // values have to be initialized
    std::cout << "\ndynamic size constructor C:\n" << dz3 << std::endl;

    // Dynamic size can be set later with resize
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>  dz4;
    dz4.resize(2,3); // values lost when resized
    dz4.setOnes(); // values have to be initialized
    std::cout << "\nDynamic resize:\n" << dz4 << std::endl;

    // resize while keeping existing values
    dz4.conservativeResize(4,3); // values not lost when resized
    std::cout << "\nDynamic resize whole keeping values:\n" << dz4 << std::endl;

    // runtime info
    std::cout << "\nruntime info:\n";
    std::cout << "size(): " << dz4.size() << std::endl;
    std::cout << "cols(): " << dz4.cols() << std::endl;
    std::cout << "rows(): " << dz4.rows() << std::endl;
    std::cout << "innerStride(): " << dz4.innerStride() << std::endl;
    std::cout << "innerSize(): " << dz4.innerSize() << std::endl;
    std::cout << "outerStride(): " << dz4.outerStride() << std::endl;
    std::cout << "outerSize(): " << dz4.outerSize() << std::endl;
    std::cout << "data(): " << dz4.data() << std::endl;

    // vector element access
    std::cout << "\nEigen::Vector4d v1(4,3,2,1)\n";
    Eigen::Vector4d v1(4,3,2,1);
    std::cout << v1(0) << " " << v1[0] << " " << v1.x() << std::endl;
    std::cout << v1.x() << " " << v1.y() << " " << v1.z() << " " << v1.w() << std::endl;
    std::cout << v1.coeff(3,0) << " " << v1(1,0) << std::endl;

    // Matrix accessing
    Eigen::MatrixXd A(4,4);
    A <<
    1, 2, 3, 4,
    5, 6, 7, 8,
    9, 8, 7, 6,
    5, 3, 2, 1;

    std::cout << "\nMatrix:\n";
    std::cout << A << std::endl;
    
    std::cout << "\nMatrix block(1, 1, 2, 2):\n";
    std::cout << A.block(1, 1, 2, 2) << std::endl;
    
    std::cout << "\nMatrix leftCols(2):\n";
    std::cout << A.leftCols(2) << std::endl;

    std::cout << "\nMatrix rightCols(2):\n";
    std::cout << A.rightCols(2) << std::endl;

    std::cout << "\nMatrix topRows(2):\n";
    std::cout << A.topRows(2) << std::endl;

    std::cout << "\nMatrix bottomRows(2):\n";
    std::cout << A.bottomRows(2) << std::endl;

    std::cout << "\nMatrix col(1):\n";
    std::cout << A.col(1) << std::endl;

    std::cout << "\nMatrix row(2):\n";
    std::cout << A.row(2) << std::endl;

    std::cout << "\nMatrix col(1).head(3):\n";
    std::cout << A.col(1).head(3) << std::endl;
    
    std::cout << "\nMatrix row(2).tail(3):\n";
    std::cout << A.row(2).tail(3) << std::endl;

    // creating column point array with dynamic count
    // transformation matrices transform these
    std::cout << "\nCreate point array:\n";
    int pointCount = 10;
    Eigen::Array4Xd  points(4,pointCount);
    points.setZero(); points.row(3).setOnes(); // make points affine

    // set value in point with index 1
    points(0,1) = 1; points(1,1) = 2; points(2,1) = 3;
    // set value in point with index 2
    // auto will create a reference of the point column, not a copy
    auto point = points.col(2);
    point(0) = 5; point(1) = 6; point(2) = 7;

    // visualization works better on transposed array
    std::cout << points.transpose() << std::endl;

    // setting point values
    std::cout << "\nChanging point values:\n";
    for ( size_t index = 0; index < points.cols(); index++ ) {
        auto point = points.col(index);
        point(0) = 6;
    }
    std::cout << points.transpose() << std::endl;

    // resize points count array while keeping values
    // points.resize(4,pointCount) will resize without keeping values
    std::cout << "\nResize point count:\n";
    pointCount = 5;
    points.conservativeResize(4,pointCount);
    std::cout << points.transpose() << std::endl;

    // Create affine transform matrix
    std::cout << "\nCreate affine transform matrix:\n";

    // set translation variable
    Eigen::Translation3d translate(2,3,4);
    
    // set rotation variables
    Eigen::AngleAxisd rotX(0, Eigen::Vector3d::UnitX());
    Eigen::AngleAxisd rotY(0, Eigen::Vector3d::UnitY());
    Eigen::AngleAxisd rotZ(M_PI, Eigen::Vector3d::UnitZ());
    
    // set Scale variables
    Eigen::AlignedScaling3d scale(2,1,1);

    // create one transformation from back to front
    // auto creates an Eigen::Transform<double,3,2>
    auto transform = translate * rotZ * rotY * rotX * scale;
    // we can turn it into a transformation matrix
    std::cout << "\nEigen::Transform<double,3,2> -> Eigen::Matrix<double, 4, 4>:\n";
    std::cout << transform.matrix().transpose() << std::endl;

    // we can also create a quaternion multiplying the AngleAxis
    // auto creates an Eigen::Quaternion<double>
    auto rotationQ = rotZ * rotY * rotX;
    std::cout << "\nEigen::Quaternion<double> from multiplying AngleAxis:\n";
    std::cout << rotationQ << std::endl;
    std::cout << "\nwEigen::Quaternion<double> -> Eigen::Matrix<double, 3, 3>:\n";
    auto rotationM = rotationQ.matrix();
    std::cout << rotationM.transpose() << std::endl;

    // Eigen::Transform can be mixed with Eigen::Quaternion
    // as long as transform is the first one it will create a Transform
    auto transformWithQ = translate * rotationQ * scale;
    // this will error
    //auto transformTestB = rotationQ * scale;
   
    // transform points with transforms
    auto transformM = transform.matrix();
    auto pointsTransform = transform.matrix() * points.matrix();
    auto transformWithQM = transformWithQ.matrix();
    auto pointsTransformWithQM = transformWithQM.matrix() * points.matrix();

    std::cout << "\nTransformed points with transform:\n";
    std::cout << pointsTransform.transpose() << std::endl;
    std::cout << "\nTransformed points with transformWithQM:\n";
    std::cout << pointsTransformWithQM.transpose() << std::endl;

    // dot and cross operations with vectors
    Eigen::Vector3d vecA = {2,0,0};
    Eigen::Vector3d vecB = {1,1,0};

    // to get the cosine value of an angle between 2 corners
    std::cout << "\nDot between 2 vectors:\n";
    auto cosVal = vecA.normalized().dot(vecB.normalized());
    auto radian = acos(cosVal);
    auto degree = ((double) 180 / M_PI) * radian;
    std::cout << "acos(" << cosVal << ") -> " << radian << " radian -> "
        << degree << " degree\n" ;

    auto crossABVal = vecA.normalized().cross(vecB.normalized());
    auto crossBAVal = vecB.normalized().cross(vecA.normalized());
    std::cout << "cross: AXB\n";
    std::cout << crossABVal.transpose() << std::endl;
    std::cout << "cross: BXA\n";
    std::cout << crossBAVal.transpose() << std::endl;
    std::cout << "sine value is in z-axis, 3rd element. Dependant of cross order.\n";
    std::cout << "resulting vector is also perpendicular to both vectors\n";

    // cross does not seem to work with dynamic sized vectors
    // have not figured out yet why
    Eigen::VectorXd vecC {{2,0,0}};
    Eigen::VectorXd vecD {{1,1,0}};
    //the following will error out
    //auto cosValX = vecC.normalized().cross(vecD.normalized());

    std::cout << "\nDone!\n";

    return 0;
}
