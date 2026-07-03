int opcion;
    long double dt_norm,t_norm;
    do
    {
        cout << "Seleccione el sistema a evaluar." << endl;
        cout << "1.- Sistema Sol-Tierra"<<endl;
        cout << "2.- Sistema Tierra-Luna"<<endl;
        cin >> opcion;
        if(opcion == 1){
             //Sistema con Sol en el origen
            Sistema.body[0].masa = aNORM_masa(1.989e30); //Sol
            Sistema.body[1].masa = aNORM_masa(5.972e24); // Tierra
            Sistema.body[2].masa = aNORM_masa(5.0e5); // Cohete 
            Sistema.body[0].posicion = {0.0,0.0,0.0};
            Sistema.body[1].posicion = {aNORM_pos(1.496e11),0,0};
            Sistema.body[2].posicion = {aNORM_pos(1.496e11+6.371e6),0.0,0.0};
            Sistema.body[0].velocidad = {0.0,0.0,0.0};
            Sistema.body[1].velocidad = {0.0,aNORM_vel(29780.0),0.0};
            Sistema.body[2].velocidad = {0.0,aNORM_vel(29780+11100.0),0.0};
    
            dt_norm = 0.01;
            t_norm = 0.0;
        }else if(opcion == 2){
            //Sistema con Tierra en el origen
            Sistema.body[0].masa = aNORM_masa(5.972e24); //Tierra
            Sistema.body[1].masa = aNORM_masa(7.348e22); // Luna
            Sistema.body[2].masa = aNORM_masa(5.0e5); // Cohete
            Sistema.body[0].posicion = {0.0,0.0,0.0};
            Sistema.body[1].posicion = {aNORM_pos(3.84e8),0.0,0.0};
            Sistema.body[2].posicion = {aNORM_pos(6.371e6),0.0,0.0};//Cohete
            Sistema.body[0].velocidad = {0.0,0.0,0.0};
            Sistema.body[1].velocidad = {0.0,aNORM_vel(1022.0),0.0};
            Sistema.body[2].velocidad = {0.0,aNORM_vel(11100.0),0.0};//Cohete
            dt_norm = 0.01;
            t_norm = 0.0;
        }else{
            cout << "Opción inválida."<<endl;
        }
    } while (opcion != 1 && opcion != 2);