#include <iostream>
#include <cmath>
#include <vector>
#include <numeric> // Para std::accumulate, std::inner_product
#include <tuple>
#include <algorithm> // Para std::transform
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

// ----------------------------------------Función para realizar la rotación de la gaussiana
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
// -------------------------- end rotacion -------------------------------------
std::pair<double, double> calcularPromedios(const std::vector<double>& x, const std::vector<double>& y) {
    // Calcular el promedio de los elementos en el vector x
    double promedio_x = std::accumulate(x.begin(), x.end(), 0.0) / x.size();
    
    // Calcular el promedio de los elementos en el vector y
    double promedio_y = std::accumulate(y.begin(), y.end(), 0.0) / y.size();
    
    return std::make_pair(promedio_x, promedio_y);
}


std::pair<std::vector<double>, std::vector<double>> dis_ang(std::vector<double> x, std::vector<double> y, double xc, double yc) {
    std::vector<double> distancias;
    std::vector<double> angulos_grados;
    
    for (size_t i = 0; i < x.size(); ++i) {
        double distancia = std::sqrt(std::pow(x[i] - xc, 2) + std::pow(y[i] - yc, 2));
        distancias.push_back(distancia);
        double angulo = std::atan2(y[i] - yc, x[i] - xc);
        double angulo_grados = std::fmod(std::fmod(std::abs(angulo) * 180 / M_PI, 360) + 360, 360);
        double equis = x[i] - xc;
        double ye = y[i] - yc;
        if (ye <= 0 && equis <= 0) {
            angulo_grados = 360 - std::abs(angulo_grados);
        } if (ye <= 0 && equis >= 0) {
            angulo_grados = 360 - std::abs(angulo_grados);
        }
        angulos_grados.push_back(angulo_grados);
    }
    
    return std::make_pair(distancias, angulos_grados);

}
//------------------------------------- FUNCION PARA ORDENAR PUNTOS -------------------------------------------------------------------------

// Función para calcular el centro de masa
std::vector<double> calcularCentroDeMasa(const std::vector<double>& x, const std::vector<double>& y) {
    double sum_x = 0, sum_y = 0;
    for (size_t i = 0; i < x.size(); ++i) {
        sum_x += x[i];
        sum_y += y[i];
    }
    return {sum_x / x.size(), sum_y / y.size()};
}

// Función para calcular los ángulos con respecto al centro de masa
std::vector<double> calcularAngulos(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& cm) {
    std::vector<double> angulos;
    for (size_t i = 0; i < x.size(); ++i) {
        angulos.push_back(std::atan2(y[i] - cm[1], x[i] - cm[0]) * 180 / M_PI);
    }
    return angulos;
}

//----------------------------------------------- END ORDENAR PUNTOS -----------------------------------------------------------------------

int main() {

    std::vector<double> x = {5000, 4000, 6000}; // Ejemplo de valores para x
    std::vector<double> y = {4000, 3000, 3000}; // Ejemplo de valores para y
//--------------------------------CONVERSION A METROS -------------------------------------------
    for (size_t i = 0; i < x.size(); ++i) {
        x[i] /= 1000; // Convertir de milímetros a metros
        y[i] /= 1000; // Convertir de milímetros a metros
    }
//---------------------------------------- ORDENAR PUNTOS ----------------------
    //Calcular el centro de masa
    std::vector<double> cm = calcularCentroDeMasa(x, y);

    // Calcular los ángulos con respecto al centro de masa (en grados)
    std::vector<double> angulos = calcularAngulos(x, y, cm);

    // Ajustar los ángulos para que estén en el rango de 0 a 360 grados
    std::transform(angulos.begin(), angulos.end(), angulos.begin(), [](double angle) {
        return std::fmod(angle + 360, 360);
    });

    // Ordenar los puntos según los ángulos ajustados
    std::vector<size_t> indices(angulos.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&](size_t i, size_t j) {
        return angulos[i] < angulos[j];
    });

    // Generar los nuevos vectores x e y ordenados
    std::vector<double> x_ordenado(x.size()), y_ordenado(y.size());
    for (size_t i = 0; i < indices.size(); ++i) {
        x_ordenado[i] = x[indices[i]];
        y_ordenado[i] = y[indices[i]];
    }

    // Mostrar los resultados
    // std::cout << "x_ordenado: ";
    // for (auto val : x_ordenado) {
    //     std::cout << val << " ";
    // }
    // std::cout << std::endl;

    // std::cout << "y_ordenado: ";
    // for (auto val : y_ordenado) {
    //     std::cout << val << " ";
    // }
    // std::cout << std::endl;

//-----------------------------------------END ORDENAR PUTNOS ------------------
    std::pair<double, double> promedios = calcularPromedios(x_ordenado, y_ordenado);
    double xc = promedios.first;
    double yc = promedios.second;
    double vf=0, vre=0, vl=0, vri=0;
    
    std::pair<std::vector<double>, std::vector<double>> resultados = dis_ang(x_ordenado, y_ordenado, xc, yc);
    std::vector<double> distancias = resultados.first;
    std::vector<double> angulos_grados = resultados.second;
    std::vector<double> sigma_x(distancias.size(), 0.0); 
    std::vector<double> sigma_y(distancias.size(), 0.0); 

    for (int i = 0; i < distancias.size(); ++i) {
        if (angulos_grados[i] >= 0 && angulos_grados[i] <= 180) {
            sigma_y[i] = std::abs((distancias[i]) * std::sin((angulos_grados[i]) * M_PI / 180)) + vf;
            sigma_x[i] = std::abs((distancias[i]) * std::cos((angulos_grados[i]) * M_PI / 180)) + vf;
        } else {
            sigma_y[i] = std::abs((distancias[i]) * std::sin((angulos_grados[i]) * M_PI / 180)) + vre;
            sigma_x[i] = std::abs((distancias[i]) * std::cos((angulos_grados[i]) * M_PI / 180)) + vre;
        }
        if (angulos_grados[i] >= 90 && angulos_grados[i] <= 270) {
            sigma_y[i] = std::abs((distancias[i]) * std::sin((angulos_grados[i]) * M_PI / 180)) + vl;
            sigma_x[i] = std::abs((distancias[i]) * std::cos((angulos_grados[i]) * M_PI / 180)) + vl;
        } else {
            sigma_y[i] = std::abs((distancias[i]) * std::sin((angulos_grados[i]) * M_PI / 180)) + vri;
            sigma_x[i] = std::abs((distancias[i]) * std::cos((angulos_grados[i]) * M_PI / 180)) + vri;
        }
    }
    sigma_x.push_back(sigma_x[0]);
    sigma_y.push_back(sigma_y[0]);

    // Determinar las distancias
    std::vector<double> distan(distancias.size(), 0.0);
    for (size_t i = 0; i < angulos_grados.size(); ++i) {
        if (i == angulos_grados.size() - 1) {
            distan[i] = 360 - angulos_grados[i] + angulos_grados[0];
        } else {
            distan[i] = angulos_grados[i + 1] - angulos_grados[i];
        }
    }

    std::vector<double> distancias_1(distan.size() + 1); // Vector para almacenar las distancias

    // Generar un valor de cero al comienzo del arreglo
    distancias_1[0] = 0;

    // Llenar el resto del arreglo con los valores de distan
    for (size_t i = 1; i < distancias_1.size(); ++i) {
        distancias_1[i] = distan[i - 1];
    }

    std::vector<double> sigma_xx(360), sigma_yy(360); // Vectores para almacenar los valores de sigma_xx y sigma_yy

    double delta_ang = 1;
    size_t j = 1;
    size_t k = 0;
    double cont = 0;
    double angulo = 0;

    for (size_t i = 0; i < 360; ++i) {
        if (i > distancias_1[j] + angulo) {
            angulo = distancias_1[j] + angulo;
            ++j;
            ++k;
            cont = 0;
        }
        double t1 = distancias_1[j] - cont;
        double t2 = cont;
        cont += delta_ang;
        sigma_xx[i] = ((t1 / distancias_1[j]) * sigma_x[k]) + ((t2 / distancias_1[j]) * sigma_x[k + 1]);
        sigma_yy[i] = ((t1 / distancias_1[j]) * sigma_y[k]) + ((t2 / distancias_1[j]) * sigma_y[k + 1]);
    }
    //hay cierta discrepancia en los valores, a estos vectores de sigma_xx y sigma_yy, sel falta el ultimo valor segun matlab.
    // generando el mesh.

    double step = 0.1; // Tamaño del paso
    double lower_limit_x = xc - 5; // Límite inferior para x
    double upper_limit_x = xc + 5; // Límite superior para x

    double lower_limit_y = yc - 5; // Límite inferior para y
    double upper_limit_y = yc + 5; // Límite superior para y

    int num_points = static_cast<int>((upper_limit_x - lower_limit_x) / step) + 1; // Número de puntos en cada dirección

    // Crear vectores para almacenar los valores de x y y
    Eigen::VectorXd x_values(num_points);
    Eigen::VectorXd y_values(num_points);

    // Llenar los vectores con valores desde lower_limit hasta upper_limit
    for (int i = 0; i < num_points; ++i) {
        x_values(i) = lower_limit_x + i * step;
        y_values(i) = lower_limit_y + i * step;
    }

    // Crear matrices meshgrid xx y yy
    Eigen::MatrixXd xx = x_values.replicate(1, num_points);
    Eigen::MatrixXd yy = y_values.transpose().replicate(num_points, 1);

    // Imprimir las dimensiones de las matrices
    std::cout << "Tamaño de la matriz xx: " << xx.rows() << "x" << xx.cols() << std::endl;
    std::cout << "Tamaño de la matriz yy: " << yy.rows() << "x" << yy.cols() << std::endl;
    yy.transposeInPlace();
    xx.transposeInPlace();

    //GRNERAR LAS VARIANZAS EN CADA PUNTO CADA 1°

    // Declarar matrices para almacenar varianzas
    Eigen::MatrixXd varianzax(xx.rows(), xx.rows());
    Eigen::MatrixXd varianzay(yy.rows(), yy.rows());

    for (int i = 0; i < xx.rows(); ++i) {
        for (int j = 0; j < xx.rows(); ++j) {
            double theta = atan2(yy(i, j) - yc, xx(i, j) - xc);
            theta = std::fmod(std::round(theta * (180 / M_PI)), 360); // Conversión a grados y módulo 360

            int alpha = static_cast<int>(theta); // Redondear al entero más cercano
            if (alpha >= 360) {
                alpha = 360;
            }
            if (alpha <= 1) {
                alpha = 1;
            }

            varianzax(i, j) = sigma_xx[alpha];
            varianzay(i, j) = sigma_yy[alpha];
        }
    }

    // Calcular la gaussiana
    Eigen::MatrixXd zz(xx.rows(), xx.cols());

    for (int i = 0; i < xx.rows(); ++i) {
        for (int j = 0; j < xx.cols(); ++j) {
            double exponente_x = std::pow(xx(i, j) - xc, 2) / (2 * std::pow(varianzax(i, j), 2));
            double exponente_y = std::pow(yy(i, j) - yc, 2) / (2 * std::pow(varianzay(i, j), 2));

            zz(i, j) = std::exp(-exponente_x - exponente_y);
        }
    }

    //determinar el conjunto de X, Y y Z que se tienen que considerar, para este caso z=0.1
    double margen = 0.1;
    double z = 0.5;
    // Crear el vector de salida
    std::vector<float> vector_x;
    std::vector<float> vector_y;
    std::vector<float> vector_z;
    // Copiar los elementos de la matriz al vector de salida
    for (int i = 0; i < xx.rows(); ++i) {
        for (int j = 0; j < xx.rows(); ++j) {

            if (zz(i,j) >= z){
                vector_z.push_back(zz(i,j));
                vector_x.push_back(xx(i,j));
                vector_y.push_back(yy(i,j));
            } 
            
        }
    }
    //----------------------------------------------------CONVERSION A MM DE LAS COORDENADAS X, Y -------------------------------
    for (size_t i = 0; i < vector_x.size(); ++i) {
        vector_x[i] *= 1000; // Convertir de milímetros a metros
        vector_y[i] *= 1000; // Convertir de milímetros a metros
    }

    //----------------------------------------END CONVERSION A MM DE LAS COORDENADAS X,Y ----------------------------------------
    std::cout << "Vector de salida: ";
    for (int i = 0; i < vector_z.size(); ++i) {
        std::cout << vector_z[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "Vector de X: ";
    for (int i = 0; i < vector_x.size(); ++i) {
        std::cout << vector_x[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "Vector de Y: ";
    for (int i = 0; i < vector_y.size(); ++i) {
        std::cout << vector_y[i] << " ";
    }
    std::cout << std::endl;
    //--------------------------------------------------------------------------------------
    // for (int i = 0; i < xx.rows(); ++i) {
    //     for (int j = 0; j < xx.cols(); ++j) {
    //             std::cout << zz(i, j) << " ";
    //         }
    //         std::cout << std::endl;
    // }
    // return 0;
}