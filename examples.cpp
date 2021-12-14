#include "Eigen/Dense"

#include <iostream>

int main()
{
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
