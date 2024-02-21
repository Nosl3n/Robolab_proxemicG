#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm> // Para std::transform
#include <numeric>   // Para std::inner_product
#include <tuple>
#include <Eigen/Dense>

// Función para convertir grados a radianes
double deg2rad(double degrees) {
    double radian = degrees * M_PI / 180.0;
    return radian;
}

// Función para determinar la traslación en los ejes
std::pair<double, double> determinarTraslacion(double cmx, double cmy) {
    double xmove;
    double ymove;
    if (cmx - cmx == 0) {
        xmove = -cmx;
    } else{
        xmove = cmx;
    }
    
    if (cmy - cmy == 0) {
        ymove = -cmy;
    } else{
        ymove = cmy;
    }
    return std::make_pair(xmove, ymove);
}

// Función para aplicar la rotación a los puntos en 3D
std::vector<std::vector<double>> aplicarRotacion(const std::vector<std::vector<double>>& puntos, double angulo) {
    // Ángulo de rotación en radianes
    angulo = deg2rad(angulo);
    int n=3;
    Eigen::MatrixXd rot (n,n);
    rot << cos(angulo), -sin(angulo), 0,
           sin(angulo),  cos(angulo), 0,
                     0,            0, 1;
    Eigen::MatrixXd punt (n,3);
    for (int i = 0; i < puntos.size(); ++i) {
        for (int j = 0; j < puntos[i].size(); ++j) {
            punt(i, j) = puntos[i][j];
        }
    }
    
    Eigen::MatrixXd puntos_rot = punt * rot;
    std::vector<std::vector<double>> rotar(punt.rows(), std::vector<double>(punt.cols()));
    for (int i = 0; i < puntos_rot.rows(); ++i) {
        for (int j = 0; j < puntos_rot.cols(); ++j) {
            rotar[i][j] = puntos_rot(i, j);
        }
    }
 
    return rotar;
}

// Función para realizar la rotación de la gaussiana
std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> rotarGaussiana(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, double angulo, double cmx, double cmy) {
    // Determinar la traslación en los ejes
    std::pair<double, double> traslacion = determinarTraslacion(cmx, cmy);
    double xmove = traslacion.first;
    double ymove = traslacion.second;
    // Aplicar la traslación de la gaussiana al origen de coordenadas
    std::vector<std::vector<double>> puntos_traslados(x.size(), std::vector<double>(3));
    for (size_t i = 0; i < x.size(); ++i) {
        puntos_traslados[i][0] = x[i] + xmove;
        puntos_traslados[i][1] = y[i] + ymove;
        puntos_traslados[i][2] = z[i];
    }
   // for (size_t i = 0; i < puntos_traslados.size(); ++i) {
     //   std::cout << "Punto " << i+1 << ": ";
       // std::cout << "X = " << puntos_traslados[i][0] << ", ";
        //std::cout << "Y = " << puntos_traslados[i][1] << ", ";
        //std::cout << "Z = " << puntos_traslados[i][2] << std::endl;
    //}
    // Aplicar la rotación a los puntos trasladados
    std::vector<std::vector<double>> puntos_rotados = aplicarRotacion(puntos_traslados, angulo);
    // Regresar la gaussiana a su posición original
    for (size_t i = 0; i < puntos_rotados.size(); ++i) {
        puntos_rotados[i][0] -= xmove;
        puntos_rotados[i][1] -= ymove;
    }
    // Separar los puntos rotados en vectores xrot, yrot, zrot
    std::vector<double> xrot(x.size());
    std::vector<double> yrot(y.size());
    std::vector<double> zrot(z.size());
    for (size_t i = 0; i < puntos_rotados.size(); ++i) {
        xrot[i] = puntos_rotados[i][0];
        yrot[i] = puntos_rotados[i][1];
        zrot[i] = puntos_rotados[i][2];
    }
    return std::make_tuple(xrot, yrot, zrot);
}

int main() {
    // Ejemplo de valores para X, Y y Z
    std::vector<double> x = {1, 2, 3};
    std::vector<double> y = {4, 5, 6};
    std::vector<double> z = {7, 8, 9};
    double angulo = 60.0;
    double cmx = 2.5;
    double cmy = 3.5;

    // Aplicar la rotación a la gaussiana
    std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> puntos_rotados = rotarGaussiana(x, y, z, angulo, cmx, cmy);

    // Obtener los puntos rotados
    std::vector<double>& xrot = std::get<0>(puntos_rotados);
    std::vector<double>& yrot = std::get<1>(puntos_rotados);
    std::vector<double>& zrot = std::get<2>(puntos_rotados);

    // Mostrar los puntos rotados
    for (size_t i = 0; i < xrot.size(); ++i) {
        std::cout << xrot[i] << " " << yrot[i] << " " << zrot[i] << std::endl;
    }

    return 0;
}
