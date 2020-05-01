//  Diego Rubi - 00121117

void createLocalA(Matrix &A,mesh m){
    float T = m.getParameter(VARIABLE_T);
    A.at(0).at(0) += -T/8;  A.at(0).at(1) += T/8;
    A.at(1).at(0) += -T/8;  A.at(1).at(1) += T/8;
}

void createLocalC(Matrix &C,mesh m){
    float ALFA = m.getParameter(VARIABLE_ALPHA);
    C.at(0).at(0) += -3.0*ALFA/2;  C.at(0).at(1) += 3.0*ALFA/2;
    C.at(1).at(0) += -3.0*ALFA/2;  C.at(1).at(1) += 3.0*ALFA/2;
}

void createLocalG(Matrix &G,mesh m){
    float LAMBDA = m.getParameter(VARIABLE_LAMBDA);
    G.at(0).at(0) += -LAMBDA/3;  G.at(0).at(1) += LAMBDA/3;
    G.at(1).at(0) += -LAMBDA/3;  G.at(1).at(1) += LAMBDA/3;
}

void createLocalI(Matrix &I,mesh m){
    float S = m.getParameter(VARIABLE_S);
    I.at(0).at(0) += -S/2.0;  I.at(0).at(1) += S/2.0;
    I.at(1).at(0) += -S/2.0;  I.at(1).at(1) += S/2.0;
}

void createLocalP(Matrix &P,mesh m){
    float l = m.getParameter(ELEMENT_LENGTH);
    float K = m.getParameter(VARIABLE_K);
    P.at(0).at(0) += K/l;      P.at(0).at(1) += -K/l;
    P.at(1).at(0) += -K/l;     P.at(1).at(1) += K/l;
}

void createLocalQ(Matrix &Q,mesh m){
    float l = m.getParameter(ELEMENT_LENGTH);
    float U = m.getParameter(VARIABLE_U);
    Q.at(0).at(0) += U/l;      Q.at(0).at(1) += -U/l;
    Q.at(1).at(0) += -U/l;     Q.at(1).at(1) += U/l;
}

Matrix createLocalK(int element,mesh &m){
    Matrix K,A,P,G,Q,C,I;

    zeroes(A,2);
    zeroes(P,2);
    zeroes(G,2);
    zeroes(Q,2);
    zeroes(C,2);
    zeroes(I,2);
    createLocalA(A,m);
    createLocalP(P,m);
    createLocalG(G,m);
    createLocalQ(Q,m);
    createLocalC(C,m);
    createLocalI(I,m);

    Vector row1, row2, row3, row4;

    row1.push_back(A.at(0).at(0)+P.at(0).at(0)); 
    row1.push_back(A.at(0).at(1)+P.at(0).at(1));
    row1.push_back(G.at(0).at(0)+Q.at(0).at(0));                  
    row1.push_back(G.at(0).at(1)+Q.at(0).at(1));

    row2.push_back(A.at(1).at(0)+P.at(1).at(0)); 
    row2.push_back(A.at(1).at(1)+P.at(1).at(1));
    row2.push_back(G.at(1).at(0)+Q.at(1).at(0)); 
    row2.push_back(G.at(1).at(1)+Q.at(1).at(1));

    row3.push_back(C.at(0).at(0)); 
    row3.push_back(C.at(0).at(1));
    row3.push_back(I.at(0).at(0)); 
    row3.push_back(I.at(0).at(1));

    row4.push_back(C.at(1).at(0)); 
    row4.push_back(C.at(1).at(1));
    row4.push_back(I.at(1).at(0)); 
    row4.push_back(I.at(1).at(1));

    K.push_back(row1); 
    K.push_back(row2); 
    K.push_back(row3); 
    K.push_back(row4);

    return K;
}

Vector createLocalb(int element,mesh &m){
    Vector b;

    float W = m.getParameter(VARIABLE_W), N = m.getParameter(VARIABLE_N), l = m.getParameter(ELEMENT_LENGTH);
    
    b.push_back(W*l/2); 
    b.push_back(W*l/2); 
    b.push_back(N*l/2); 
    b.push_back(N*l/2);

    return b;
}

void crearSistemasLocales(mesh &m,vector<Matrix> &localKs,vector<Vector> &localbs){
    for(int i=0;i<m.getSize(ELEMENTS);i++){
        localKs.push_back(createLocalK(i,m));
        localbs.push_back(createLocalb(i,m));
    }
}

void assemblyK(element e,Matrix localK,Matrix &K,int nnodes){
    int index1 = e.getNode1() - 1;
    int index2 = e.getNode2() - 1;
    int index3 = index1 + nnodes;
    int index4 = index2 + nnodes;

    K.at(index1).at(index1) += localK.at(0).at(0);
    K.at(index1).at(index2) += localK.at(0).at(1);
    K.at(index2).at(index1) += localK.at(1).at(0);
    K.at(index2).at(index2) += localK.at(1).at(1);

    K.at(index1).at(index3) += localK.at(0).at(2);
    K.at(index1).at(index4) += localK.at(0).at(3);
    K.at(index2).at(index3) += localK.at(1).at(2);
    K.at(index2).at(index4) += localK.at(1).at(3);

    K.at(index3).at(index1) += localK.at(2).at(0);
    K.at(index3).at(index2) += localK.at(2).at(1);
    K.at(index4).at(index1) += localK.at(3).at(0);
    K.at(index4).at(index2) += localK.at(3).at(1);

    K.at(index3).at(index3) += localK.at(2).at(2);
    K.at(index3).at(index4) += localK.at(2).at(3);
    K.at(index4).at(index3) += localK.at(3).at(2);
    K.at(index4).at(index4) += localK.at(3).at(3);

}

void assemblyb(element e,Vector localb,Vector &b){
    int index1 = e.getNode1() - 1;
    int index2 = e.getNode2() - 1;

    b.at(index1) += localb.at(0);
    b.at(index2) += localb.at(1);
}

void ensamblaje(mesh &m,vector<Matrix> &localKs,vector<Vector> &localbs,Matrix &K,Vector &b){
    int nnodes = m.getSize(NODES);
    for(int i=0;i<m.getSize(ELEMENTS);i++){
        element e = m.getElement(i);
        assemblyK(e,localKs.at(i),K,nnodes);
        assemblyb(e,localbs.at(i),b);
    }
}


void applyDirichlet(mesh &m,Matrix &K,Vector &b){
    for(int i=0;i<m.getSize(DIRICHLET);i++){
        condition c = m.getCondition(i,DIRICHLET);

        int index = c.getNode1()-1;

        K.erase(K.begin()+index);
        b.erase(b.begin()+index);

        for(int row=0;row<K.size();row++){
            float cell = K.at(row).at(index);
            K.at(row).erase(K.at(row).begin()+index);
            b.at(row) += -1*c.getValue()*cell;
        }
    }
}


void calculate(Matrix &K, Vector &b, Vector &T){
    cout << "Iniciando calculo de respuesta...\n";
    Matrix Kinv;
    cout << "Calculo de inversa...\n";
    inverseMatrix(K,Kinv);
    cout << "Calculo de respuesta...\n";
    productMatrixVector(Kinv,b,T);
}
