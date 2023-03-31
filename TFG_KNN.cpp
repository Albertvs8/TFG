#define PROFILE

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "openfhe.h"
#include <cmath>

using namespace lbcrypto;
using namespace std;


//PAS 4
Ciphertext<DCRTPoly> f_4_homomorphic(const CryptoContext<DCRTPoly> &cc, const Ciphertext<DCRTPoly> &c, int vec_size) {
    auto c2 = cc->EvalMult(c, c); // c^2
    auto c3 = cc->EvalMult(c2, c); // c^3
    auto c5 = cc->EvalMult(c2, c3); // c^5
    auto c7 = cc->EvalMult(c2, c5); // c^7
    auto c9 = cc->EvalMult(c2, c7); // c^9

    auto const1 = cc->MakeCKKSPackedPlaintext(std::vector<double>(vec_size,{35.0/128.0}));
    auto const2 = cc->MakeCKKSPackedPlaintext(std::vector<double>(vec_size,{-180.0/128.0}));
    auto const3 = cc->MakeCKKSPackedPlaintext(std::vector<double>(vec_size,{378.0/128.0}));
    auto const4 = cc->MakeCKKSPackedPlaintext(std::vector<double>(vec_size,{-420.0/128.0}));
    auto const5 = cc->MakeCKKSPackedPlaintext(std::vector<double>(vec_size,{315.0/128.0}));


    auto term1 = cc->EvalMult(const1, c9);
    auto term2 = cc->EvalMult(const2, c7);
    auto term3 = cc->EvalMult(const3, c5);
    auto term4 = cc->EvalMult(const4, c3);
    auto term5 = cc->EvalMult(const5, c);
    auto result = cc->EvalAdd(term1, term2);
    result = cc->EvalAdd(result, term3);
    result = cc->EvalAdd(result, term4);
    result = cc->EvalAdd(result, term5);

    return result;
}

Ciphertext<DCRTPoly> g_4_homomorphic(const CryptoContext<DCRTPoly> &cc, const Ciphertext<DCRTPoly> &c, int vec_size) {
    auto c2 = cc->EvalMult(c, c); // c^2
    auto c3 = cc->EvalMult(c2, c); // c^3
    auto c5 = cc->EvalMult(c2, c3); // c^5
    auto c7 = cc->EvalMult(c2, c5); // c^7
    auto c9 = cc->EvalMult(c2, c7); // c^9

    auto const1 = cc->MakeCKKSPackedPlaintext(std::vector<double>(vec_size,{46623.0 / 1024.0}));
    auto const2 = cc->MakeCKKSPackedPlaintext(std::vector<double>(vec_size,{-113492.0 / 1024.0}));
    auto const3 = cc->MakeCKKSPackedPlaintext(std::vector<double>(vec_size,{97015.0 / 1024.0}));
    auto const4 = cc->MakeCKKSPackedPlaintext(std::vector<double>(vec_size,{-34974.0 / 1024.0}));
    auto const5 = cc->MakeCKKSPackedPlaintext(std::vector<double>(vec_size,{5850.0 / 1024.0}));


    auto term1 = cc->EvalMult(const1, c9);
    auto term2 = cc->EvalMult(const2, c7);
    auto term3 = cc->EvalMult(const3, c5);
    auto term4 = cc->EvalMult(const4, c3);
    auto term5 = cc->EvalMult(const5, c);
    auto result = cc->EvalAdd(term1, term2);
    result = cc->EvalAdd(result, term3);
    result = cc->EvalAdd(result, term4);
    result = cc->EvalAdd(result, term5);

    return result;
}



Ciphertext<DCRTPoly> homomorphic_comparison_g(const CryptoContext<DCRTPoly> &cc, const Ciphertext<DCRTPoly> &c1, const Ciphertext<DCRTPoly> &c2, int df, int dg, int vec_size) {
    auto x = cc->EvalSub(c1, c2);
    for (int i=1; i<=dg; i++){
        x = g_4_homomorphic(cc,x,vec_size);
    }
    for (int i=1; i<=df; i++){
        x = f_4_homomorphic(cc,x,vec_size);
    }
    auto homomorphic_1 = cc->MakeCKKSPackedPlaintext(std::vector<double>(vec_size,{1.0}));
    auto homomorphic_05 = cc->MakeCKKSPackedPlaintext(std::vector<double>(vec_size,{0.5}));
    auto x_plus_1 = cc->EvalAdd(x, homomorphic_1);
    auto result = cc->EvalMult(x_plus_1, homomorphic_05);
    return result;
}

vector<Ciphertext<DCRTPoly>> L_function(const CryptoContext<DCRTPoly> &cc, const Ciphertext<DCRTPoly> comparison, const vector<Ciphertext<DCRTPoly>> &F, const vector<Ciphertext<DCRTPoly>> &G ){
    //L_(a>b)(F,G) = (a>b)*F + (a<b)*G
 
    //L function calculation for distance
    auto first_term_dist = cc->EvalMult(comparison,F[0]);
    auto second_term_dist = cc->EvalMult(cc->EvalAdd(1,-comparison),G[0]);
    auto result_dist = cc->EvalAdd(first_term_dist,second_term_dist);

    //L function calculation for label
    auto first_term_label = cc->EvalMult(comparison,F[1]);
    auto second_term_label = cc->EvalMult(cc->EvalAdd(1,-comparison),G[1]);
    auto result_label = cc->EvalAdd(first_term_label,second_term_label);

    vector<Ciphertext<DCRTPoly>> result = {result_dist,result_label};
    return result;

}

vector<Ciphertext<DCRTPoly>> m_max(int m, const CryptoContext<DCRTPoly> &cc, const vector<vector<Ciphertext<DCRTPoly>>> &B, const vector<vector<Ciphertext<DCRTPoly>>> &C,const vector<vector<Ciphertext<DCRTPoly>>> &comparisonSetBC) 
{   
    int size_b = B.size();
    int size_c = C.size();

    //base case
    if (size_b == 0 || size_c == 0) {
       return size_b == 0 ? C[m - 1] : B[m - 1];
    } 

    //base case, when M=1 return max of first elements of both arrays
    if (m==1) {
        return L_function(cc, comparisonSetBC[0][0], B[0], C[0]);
    } 
    
    //base case, when M=size_b+size_c return min of last elements of both arrays
    if (m==size_b+size_c){
        return L_function(cc, comparisonSetBC[size_b-1][size_c-1], C[size_c-1], B[size_b-1]);
    }


    int i = floor(m / 2.0);
    int j = ceil(m / 2.0);

    vector<Ciphertext<DCRTPoly>> left;
    vector<Ciphertext<DCRTPoly>> right;

    int n=comparisonSetBC.size();
    int p=comparisonSetBC[0].size();

    if (i<size_b){ // i<=size_b
        if (j<=size_c)
        {
            vector<vector<Ciphertext<DCRTPoly>>> comparisonSetLEFT(n-i);
            for (int x = 0; x < n-i; x++) {
                comparisonSetLEFT[x].resize(j);
                for (int y = 0; y < j; y++) {
                    comparisonSetLEFT[x][y] = comparisonSetBC[i+x][y];
                }
            }
            left = m_max(j, cc, vector<vector<Ciphertext<DCRTPoly>>>(B.begin() + i, B.end()), vector<vector<Ciphertext<DCRTPoly>>>(C.begin(), C.begin() + j),comparisonSetLEFT);

            vector<vector<Ciphertext<DCRTPoly>>> comparisonSetRIGHT(i);
            for (int x = 0; x < i; x++) {
                comparisonSetRIGHT[x].resize(p - j);
                for (int y = 0; y < p-j; y++) {
                    comparisonSetRIGHT[x][y] = comparisonSetBC[x][j+y];
                }
            }
            right = m_max(i, cc, vector<vector<Ciphertext<DCRTPoly>>>(B.begin(), B.begin() + i), vector<vector<Ciphertext<DCRTPoly>>>(C.begin() + j, C.end()),comparisonSetRIGHT);
        } 
        
        else 
        {
            vector<vector<Ciphertext<DCRTPoly>>> comparisonSetLEFT(n-i);
            for (int x = 0; x < n-i; x++) {
                comparisonSetLEFT[x].resize(p);
                for (int y = 0; y < p; y++) {
                    comparisonSetLEFT[x][y] = comparisonSetBC[i+x][y];
                }
            }
            left = m_max(j, cc, vector<vector<Ciphertext<DCRTPoly>>>(B.begin() + i, B.end()), vector<vector<Ciphertext<DCRTPoly>>>(C.begin(), C.end()),comparisonSetLEFT);

            vector<vector<Ciphertext<DCRTPoly>>> comparisonSetRIGHT = vector<vector<Ciphertext<DCRTPoly>>>(comparisonSetBC.begin(), comparisonSetBC.begin());
            right = m_max(i, cc, vector<vector<Ciphertext<DCRTPoly>>>(B.begin(), B.begin() + i), vector<vector<Ciphertext<DCRTPoly>>>(C.begin(), C.begin()),comparisonSetRIGHT);
        } 
    }
    else {

        vector<vector<Ciphertext<DCRTPoly>>> comparisonSetLEFT = vector<vector<Ciphertext<DCRTPoly>>>(comparisonSetBC.begin(), comparisonSetBC.begin());
        left = m_max(j, cc, vector<vector<Ciphertext<DCRTPoly>>>(B.begin(), B.begin()), vector<vector<Ciphertext<DCRTPoly>>>(C.begin(), C.begin() + j),comparisonSetLEFT);

        vector<vector<Ciphertext<DCRTPoly>>> comparisonSetRIGHT(n);
            for (int x = 0; x < n; x++) {
                comparisonSetRIGHT[x].resize(p - j);
                for (int y = 0; y < p-j; y++) {
                    comparisonSetRIGHT[x][y] = comparisonSetBC[x][j+y];
                }
            }
        right = m_max(i, cc, vector<vector<Ciphertext<DCRTPoly>>>(B.begin(), B.end()), vector<vector<Ciphertext<DCRTPoly>>>(C.begin() + j, C.end()),comparisonSetRIGHT);
    }

    Ciphertext<DCRTPoly> xiCOMPyj;
    if (i>size_b){
        xiCOMPyj = comparisonSetBC[size_b-1][j-1];
    } else {
        if (j>size_c){
            xiCOMPyj = comparisonSetBC[i-1][size_c-1];
        }
        else{
            xiCOMPyj = comparisonSetBC[i-1][j-1];
        }
    }
    return L_function(cc, xiCOMPyj, left, right);
}



vector<vector<Ciphertext<DCRTPoly>>> merge(const CryptoContext<DCRTPoly> &cc, const vector<vector<Ciphertext<DCRTPoly>>> &B, const vector<vector<Ciphertext<DCRTPoly>>> &C, const vector<vector<Ciphertext<DCRTPoly>>> &comparisonSetBC){
    int s = B.size();
    int t = C.size();
    //int k = floor((s+t)/2);


    vector<vector<Ciphertext<DCRTPoly>>> Z(s+t);
    for (int i=1; i<=(s+t); i++){
        vector<Ciphertext<DCRTPoly>> zi = m_max(i,cc,B,C,comparisonSetBC);
        Z[i-1] = zi;
        
    }
    return Z;
}


vector<vector<Ciphertext<DCRTPoly>>> merge_sort(const CryptoContext<DCRTPoly> &cc, const vector<vector<Ciphertext<DCRTPoly>>> &A, const vector<vector<Ciphertext<DCRTPoly>>> &comparisonSetA, const Ciphertext<DCRTPoly> &aux_enc1){
    if (A.size()==1) {
        return A;
    }
    int s=floor(A.size()/2);
    int n=A.size(); 
    

    vector<vector<Ciphertext<DCRTPoly>>> comparisonSetB(s);

    for (int i = 0; i < s; i++) {
        comparisonSetB[i].resize(s);
        for (int j = 0; j < s; j++) {
            comparisonSetB[i][j] = comparisonSetA[i][j];
        }
    }

    vector<vector<Ciphertext<DCRTPoly>>> comparisonSetC(n-s); 
    for (int i = 0; i < n-s; i++) {
        comparisonSetC[i].resize(n-s);
        for (int j = 0; j < n-s; j++) {
            comparisonSetC[i][j] = comparisonSetA[s+i][s+j];
        }
    }

      
    auto B = merge_sort(cc,vector<vector<Ciphertext<DCRTPoly>>>(A.begin(), A.begin()+s), comparisonSetB,aux_enc1);
    auto C = merge_sort(cc,vector<vector<Ciphertext<DCRTPoly>>>(A.begin()+s, A.end()), comparisonSetC,aux_enc1);

    int k=A.size();

    vector<vector<vector<Ciphertext<DCRTPoly>>>> unionBiAj;
    for (int j = (s+1); j<=k; j++){
        auto sorterParameter = vector<vector<Ciphertext<DCRTPoly>>>(comparisonSetA.begin(), comparisonSetA.begin()+s); //ERROR AQUI
        vector<vector<Ciphertext<DCRTPoly>>> sorterParameter1;
        int sp_size = sorterParameter.size();

        for (int i=0; i<sp_size; i++){
            vector<Ciphertext<DCRTPoly>> row = {sorterParameter[i][j-1], aux_enc1};
            sorterParameter1.push_back(row);
        }

        auto BCOMPaj = merge_sort(cc,sorterParameter1,comparisonSetB,aux_enc1);
        unionBiAj.push_back(BCOMPaj);
    }



    vector<vector<vector<Ciphertext<DCRTPoly>>>> unionBC;
    for (int i = 1; i<=s; i++){
        vector<vector<Ciphertext<DCRTPoly>>> sorterParameter2;
        for(int j=0; j<k-s; j++){ 
            sorterParameter2.push_back(unionBiAj[j][i-1]);
        }
        auto biCOMPC = merge_sort(cc,sorterParameter2,comparisonSetC,aux_enc1);
        unionBC.push_back(biCOMPC);
    }




    vector<vector<Ciphertext<DCRTPoly>>> BCompC;
    int size_unionBC = unionBC.size();
    int size_inner_unionBC = unionBC[0].size();


    BCompC.resize(size_unionBC);
    for (int i=0; i<size_unionBC; i++){
        BCompC[i].resize(size_inner_unionBC);
        for (int j=1; j<=size_inner_unionBC; j++){
            BCompC[i][j-1] = unionBC[i][j-1][0];
        }
    }

    auto return_value = merge(cc,B,C,BCompC);
    return return_value;

}

int main() {
    
    
    //STEP 0 - Llegir dades
     // Create a vector to store the dataset
    vector<vector<double>> dataset;

    vector<double> row1 = {25.5, 1.6, 0.8, 0}; // Patient age, BMI, blood sugar level, has diabetes (0=no, 1=yes)
    vector<double> row2 = {32.0, 2.2, 1.3, 1};
    vector<double> row3 = {45.3, 3.1, 1.7, 1};
    vector<double> row4 = {38.9, 2.8, 1.4, 1};
    vector<double> row5 = {29.1, 1.9, 0.9, 0};
    vector<double> row6 = {50.2, 3.6, 1.9, 1};

    dataset.push_back(row1);
    dataset.push_back(row2);
    dataset.push_back(row3);
    dataset.push_back(row4);
    dataset.push_back(row5);
    dataset.push_back(row6);

    // Print the dataset with headers
    cout << "Age\tBMI\tBS\tDiabetes" << endl;
    for (size_t i = 0; i < dataset.size(); i++) {
        for (size_t j = 0; j < dataset[i].size(); j++) {
            cout << dataset[i][j] << "\t";
        }
        cout << endl;
    }




    //STEP 1 - Escalar perque les distàncies estiguin entre 0 i 1


    //STEP 2- Encriptar dades



    //STEP 3 - Calcular euclidean distance (per cada sample de test)

    //STEP 4 - Ordenar euclidean distance (per cada sample de test), fent servir merge_sort
    /*int df=2;
    int dg=2;
    int vec_size=1;*/


    //STEP 5 - Agafar labels de les k distàncies més petites

    //STEP 6 - Sumar les k labels


    //STEP 7 - Desencriptar resultat de step 6) i en funció d'això decidir la label

    return 0;
}




    
    


    

    
   

    
