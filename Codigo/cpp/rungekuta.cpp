#include <bits/stdc++.h>
using namespace std;
//Constante de gravitacion universal
const double G = 6.67430e-11;  // m^3 kg^-1 s^-2
//-----------------------------------------------------
// Definimos nuestros tres cuerpos y vectores tridimensionales
//-----------------------------------------------------
const int N = 3; 
struct Vector3D
{
    double x, y,z;

    Vector3D operator+(const Vector3D& v) const
    {
        return {x + v.x, y + v.y,z+v.z};
    }

    Vector3D operator-(const Vector3D& v) const
    {
        return {x - v.x, y - v.y,z-v.z};
    }

    Vector3D operator*(double k) const
    {
        return {k*x, k*y,k*z};
    }
     Vector3D operator/(double k) const
    {
        return {x/k, y/k, z/k};
    }

    double norma() const
    {
        return sqrtl(x*x + y*y + z*z);
    }

    double normaAlCuadrado() const
    {
        return x*x + y*y + z*z;
    }

    Vector3D normalizar() const
    {
        return (*this)/norma();
    }
    Vector3D& operator+=(const Vector3D& v)
    {
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }
};
struct Body{
    double masa;
    Vector3D posicion;
    Vector3D velocidad;
};
struct State
{
    Body body[N];
};
struct Derivative
{
    Vector3D dPosition[N];

    Vector3D dVelocity[N];
};
double productoInterno(Vector3D v , Vector3D w){
    return v.x*w.x+v.y*w.y+v.z*w.z;
}
Vector3D productoCruz(Vector3D v , Vector3D w ){
    return {v.y*w.z-w.y*v.z,v.z*w.x-v.x*w.z,v.x*w.y-w.x*v.y};
}

Derivative evaluate(const State& state)
{
    Derivative k;

    // Las derivadas de la posición son las velocidades
    for(int i=0;i<N;i++)
        k.dPosition[i] = state.body[i].velocidad;

    // Inicializar aceleraciones
    for(int i=0;i<N;i++)
        k.dVelocity[i] = {0.0,0.0,0.0};

    // Calcular aceleraciones
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            if(i==j){
                continue;
            }
            Vector3D r =
                state.body[j].posicion
              - state.body[i].posicion;

            double d = r.norma();

            double factor =
                G*state.body[j].masa/(d*d*d);

            k.dVelocity[i] += r*factor;
        }
    }

    return k;
}

State addState(const State& state,
               const Derivative& k,
               double h)
{
    State result = state;

    for(int i=0;i<N;i++)
    {
        result.body[i].posicion +=
            k.dPosition[i]*h;

        result.body[i].velocidad +=
            k.dVelocity[i]*h;
    }

    return result;
}
//Metodo Runge Kuta de 4to orden
void RK4Step(State& state,double h)
{
    Derivative k1 = evaluate(state);

    Derivative k2 =
        evaluate(addState(state,k1,h/2.0));

    Derivative k3 =
        evaluate(addState(state,k2,h/2.0));

    Derivative k4 =
        evaluate(addState(state,k3,h));

    for(int i=0;i<N;i++)
    {
        state.body[i].posicion +=
        (k1.dPosition[i]
        +k2.dPosition[i]*2.0
        +k3.dPosition[i]*2.0
        +k4.dPosition[i])*(h/6.0);

        state.body[i].velocidad +=
        (k1.dVelocity[i]
        +k2.dVelocity[i]*2.0
        +k3.dVelocity[i]*2.0
        +k4.dVelocity[i])*(h/6.0);
    }
}

double Energia(State& state){
    double kinetic = 0.0;
    double potential = 0.0;
    Vector3D temp;
    for (int i = 0; i < N; i++)
    {
        kinetic+= state.body[i].masa * state.body[i].velocidad.normaAlCuadrado()/2;
        for (int j = i+1; j < N; j++)
        {
            temp = state.body[i].posicion-state.body[j].posicion;
            potential-=G*state.body[i].masa * state.body[j].masa/temp.norma();
        }
    }
    return potential + kinetic;
}
Vector3D momentoLineal(State& state){
    Vector3D momentoResultante = {0.0,0.0,0.0};
    for (int i = 0; i < N; i++)
    {
        momentoResultante += state.body[i].velocidad*state.body[i].masa;
    }
    return momentoResultante;
}

Vector3D momentoAngular(State& state){
    Vector3D momentoAngularResultante = {0.0,0.0,0.0};
    for (int i = 0; i < N; i++)
    {
        momentoAngularResultante+= productoCruz(state.body[i].posicion,state.body[i].velocidad)*state.body[i].masa;
    }
    return momentoAngularResultante;
}
//-----------------------------------------------------
// Programa principal
//-----------------------------------------------------
State Sistema;
int main()
{
    // Inicializar masas
    fstream archivo;
    fstream archivo2;
    fstream archivo3;
    fstream archivo4;
    archivo.open("rungekuta.csv",ios::in);
    if(archivo.good()){
        archivo.close();
        cout << "El archivo se abrió con éxito"<<endl;
    }
    int opcion;
    double dt,t;
    do
    {
        cout << "Seleccione el sistema a evaluar." << endl;
        cout << "1.- Sistema Sol-Tierra"<<endl;
        cout << "2.- Sistema Tierra-Luna"<<endl;
        cin >> opcion;
        if(opcion == 1){
             //Sistema con Sol en el origen
            Sistema.body[0].masa = 1.989e30; //Sol
            Sistema.body[1].masa = 5.972e24; // Tierra
            Sistema.body[2].masa = 5.0e5; // Cohete 
            Sistema.body[0].posicion = {0.0,0.0,0.0};
            Sistema.body[1].posicion = {1.496e11,0,0};
            Sistema.body[2].posicion = {1.496e11+6.371e6,0.0,0.0};
            Sistema.body[0].velocidad = {0.0,0.0,0.0};
            Sistema.body[1].velocidad = {0.0,29780.0,0.0};
            Sistema.body[2].velocidad = {0.0,29780+11100.0,0.0};
    
            dt = 80.5469;;
            t = 0.0;
        }else if(opcion == 2){
            //Sistema con Tierra en el origen
            Sistema.body[0].masa = 5.972e24; //Tierra
            Sistema.body[1].masa = 7.348e22; // Luna
            Sistema.body[2].masa = 5.0e5; // Cohete
            Sistema.body[0].posicion = {0.0,0.0,0.0};
            Sistema.body[1].posicion = {3.84e8,0.0,0.0};
            Sistema.body[2].posicion = {6.371e6,0.0,0.0};//Cohete
            Sistema.body[0].velocidad = {0.0,0.0,0.0};
            Sistema.body[1].velocidad = {0.0,1022.0,0.0};
            Sistema.body[2].velocidad = {0.0,11100.0,0.0};//Cohete
            dt = 8.05469;;
            t = 0.0;
        }else{
            cout << "Opción inválida."<<endl;
        }
    } while (opcion != 1 && opcion != 2);
    //Sistema con SOl en el origen
    /* Sistema.body[0].masa = 1.989e30; //Sol
    Sistema.body[1].masa = 5.972e24; // Tierra
    //Sistema.body[2].masa = 7.348e22; // Luna
    Sistema.body[2].masa = 5.0e5; // Cohete */
    /* Sistema.body[0].masa = 1.989e30; //Sol
     Sistema.body[1].masa = 5.972e24; // Tierra
     Sistema.body[2].masa = 5.0e5; // Cohete */
  

    // Inicializar posiciones
    //Sistema con Sol en el origen
   /*  Sistema.body[0].posicion = {0.0,0.0,0.0};
     Sistema.body[1].posicion = {1.496e11,0,0};
    Sistema.body[2].posicion = {1.496e11+6.371e6,0.0,0.0}; */
    // Sistema con Tierra en el origen

    //Sistema.body[2].posicion = {1.496e11+3.84e8,0,0};//Luna

   /*  Sistema.body[0].velocidad = {0.0,0.0,0.0};
    Sistema.body[1].velocidad = {0.0,29780.0,0.0};
    Sistema.body[2].velocidad = {0.0,29780+11100.0,0.0}; */
    //Sistema con Tierra en el origen

    //Sistema.body[2].velocidad = {0,29780+1022,0};//Luna
    const int cantidadDePasos = 1000000;
     const int GUARDAR_CADA = 100;
    archivo.open("rungekuta.csv", ios::out);
    archivo2.open("sistemalunatierra.csv", ios::out);
    archivo3.open("sistemasateliteluna.csv", ios::out);
    archivo4.open("cantidadesconservadas.csv", ios::out);
    archivo << "t,"
     << "x1,y1,z1,"
     << "x2,y2,z2,"
     << "x3,y3,z3\n";
    archivo2 << "t,"
     << "x1,y1,z1\n";
    archivo3 << "t,"
     << "x1,y1,z1\n";
    archivo4 << "t,"
    <<"E,"
    <<"p1,p2,p3,"
    <<"L1,L2,L3\n";
    for(int i=0;i<cantidadDePasos;i++)
    {
        RK4Step(Sistema,dt);
        //FORMATO
        if(i%GUARDAR_CADA == 0){
            Vector3D mlineal = momentoLineal(Sistema);
            Vector3D mangular = momentoAngular(Sistema);
        archivo << t <<","
        << Sistema.body[0].posicion.x << ","
        << Sistema.body[0].posicion.y << ","
        << Sistema.body[0].posicion.z << ","
        << Sistema.body[1].posicion.x << ","
        << Sistema.body[1].posicion.y << ","
        << Sistema.body[1].posicion.z << ","
        << Sistema.body[2].posicion.x << ","
        << Sistema.body[2].posicion.y << ","
        << Sistema.body[2].posicion.z << "\n";
        archivo2 <<t << ","
        << Sistema.body[2].posicion.x-Sistema.body[0].posicion.x << ","
        << Sistema.body[2].posicion.y-Sistema.body[0].posicion.y << ","
        << Sistema.body[2].posicion.z-Sistema.body[0].posicion.z << "\n";
        archivo3 <<t << ","
        << Sistema.body[2].posicion.x-Sistema.body[1].posicion.x << ","
        << Sistema.body[2].posicion.y-Sistema.body[1].posicion.y << ","
        << Sistema.body[2].posicion.z-Sistema.body[1].posicion.z << "\n";
        archivo4 << t << ","
        << Energia(Sistema) <<","
        <<mlineal.x << ","
        << mlineal.y << ","
        << mlineal.z << ","
        << mangular.x << ","
        << mangular.y << ","
        << mangular.z << "\n";
        }
        t+=dt;
    }
    
    archivo.close();
    archivo2.close();
    archivo3.close();
    archivo4.close();
    return 0;
}