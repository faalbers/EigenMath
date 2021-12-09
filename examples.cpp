#include "Eigen/Dense"

#include <iostream>

int main()
{
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

    // creating affine transformation
    Eigen::Transform<double,3,Eigen::Affine> transform;
    // transformations go from back to front
    transform = translate * rotZ * rotY * rotX * scale;
    std::cout << transform.matrix().transpose() << std::endl;

    // transform points with transform 
    points = transform.matrix() * points.matrix();

    std::cout << "\nTransformed points:\n";
    std::cout << points.transpose() << std::endl;

    std::cout << "\nDone!\n";

    return 0;
}
