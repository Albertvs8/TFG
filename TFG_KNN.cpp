#define PROFILE

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "openfhe.h"
#include <cmath>
#include <stdexcept>

using namespace lbcrypto;
using namespace std;


//PAS 1
// Function to scale the dataset using min-max scaling
vector<vector<vector<double>>> min_max_scaling(const vector<vector<vector<double>>>& dataset) {
    vector<vector<double>> min_max_values(dataset[0].size()-1, vector<double>(2));
    for (auto& row : min_max_values) {
        row[0] = std::numeric_limits<double>::max();
        row[1] = std::numeric_limits<double>::lowest();
    }

    for (const auto& row : dataset) {
        for (size_t i = 0; i < row.size()-1; i++) {
            double value = row[i][0];
            if (value < min_max_values[i][0]) {
                min_max_values[i][0] = value; 
            }
            if (value > min_max_values[i][1]) {
                min_max_values[i][1] = value; 
            }
        }
    }
    
    vector<vector<vector<double>>> scaled_dataset;
    for (const auto& row : dataset) {
        vector<vector<double>> scaled_row(row.size(), vector<double>(1));
        for (size_t i = 0; i < row.size()-1; i++) {
            double value = row[i][0];
            double min_value = min_max_values[i][0];
            double max_value = min_max_values[i][1];
            double scaled_value = (0.57*(value - min_value)) / (max_value - min_value);
            scaled_row[i][0] = scaled_value;
        }
        scaled_row.back()[0] = row.back()[0];
        scaled_dataset.push_back(scaled_row);
    }

    return scaled_dataset;
}



//PAS 3
Ciphertext<DCRTPoly> euclidean_distance(const CryptoContext<DCRTPoly> &cc, const vector<Ciphertext<DCRTPoly>> &X, const vector<Ciphertext<DCRTPoly>> &Y){
    if (X.size()!=Y.size()){
        throw std::invalid_argument("X and Y sizes are different");
    }

    vector<Ciphertext<DCRTPoly>> potencies;
    for (size_t i = 0; i < X.size(); i++) { 
        auto resta = cc->EvalAdd(X[i],-Y[i]);
        auto potencia = cc->EvalMult(resta, resta);
        potencies.push_back(potencia);
    }

    Ciphertext<DCRTPoly> result = potencies[0];
    for (size_t i = 1; i < X.size(); i++) { 
        result = cc->EvalAdd(result,potencies[i]);
    }
    return result;
}


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

Ciphertext<DCRTPoly> f_3_homomorphic(const CryptoContext<DCRTPoly> &cc, const Ciphertext<DCRTPoly> &c, int vec_size) {
    auto c2 = cc->EvalMult(c, c); // c^2
    auto c3 = cc->EvalMult(c2, c); // c^3
    auto c5 = cc->EvalMult(c2, c3); // c^5
    auto c7 = cc->EvalMult(c2, c5); // c^7

    auto const1 = cc->MakeCKKSPackedPlaintext(std::vector<double>(vec_size,{-5.0/16.0}));
    auto const2 = cc->MakeCKKSPackedPlaintext(std::vector<double>(vec_size,{21.0/16.0}));
    auto const3 = cc->MakeCKKSPackedPlaintext(std::vector<double>(vec_size,{-35.0/16.0}));
    auto const4 = cc->MakeCKKSPackedPlaintext(std::vector<double>(vec_size,{35.0/16.0}));


    auto term1 = cc->EvalMult(const1, c7);
    auto term2 = cc->EvalMult(const2, c5);
    auto term3 = cc->EvalMult(const3, c3);
    auto term4 = cc->EvalMult(const4, c);

    auto result = cc->EvalAdd(term1, term2);
    result = cc->EvalAdd(result, term3);
    result = cc->EvalAdd(result, term4);

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

Ciphertext<DCRTPoly> g_3_homomorphic(const CryptoContext<DCRTPoly> &cc, const Ciphertext<DCRTPoly> &c, int vec_size) {
    auto c2 = cc->EvalMult(c, c); // c^2
    auto c3 = cc->EvalMult(c2, c); // c^3
    auto c5 = cc->EvalMult(c2, c3); // c^5
    auto c7 = cc->EvalMult(c2, c5); // c^7

    auto const1 = cc->MakeCKKSPackedPlaintext(std::vector<double>(vec_size,{-12860.0/1024.0}));
    auto const2 = cc->MakeCKKSPackedPlaintext(std::vector<double>(vec_size,{25614.0/1024.0}));
    auto const3 = cc->MakeCKKSPackedPlaintext(std::vector<double>(vec_size,{-16577.0/1024.0}));
    auto const4 = cc->MakeCKKSPackedPlaintext(std::vector<double>(vec_size,{4589.0/1024.0}));


    auto term1 = cc->EvalMult(const1, c7);
    auto term2 = cc->EvalMult(const2, c5);
    auto term3 = cc->EvalMult(const3, c3);
    auto term4 = cc->EvalMult(const4, c);

    auto result = cc->EvalAdd(term1, term2);
    result = cc->EvalAdd(result, term3);
    result = cc->EvalAdd(result, term4);

    return result;
}



Ciphertext<DCRTPoly> homomorphic_comparison_g(const CryptoContext<DCRTPoly> &cc, const Ciphertext<DCRTPoly> &c1, const Ciphertext<DCRTPoly> &c2, int df, int dg, int vec_size) {
    auto x = cc->EvalSub(c1, c2);
    for (int i=1; i<=dg; i++){
        x = g_3_homomorphic(cc,x,vec_size);
    }
    for (int i=1; i<=df; i++){
        x = f_3_homomorphic(cc,x,vec_size);
    }
    auto homomorphic_1 = cc->MakeCKKSPackedPlaintext(std::vector<double>(vec_size,{1.0}));
    auto homomorphic_05 = cc->MakeCKKSPackedPlaintext(std::vector<double>(vec_size,{0.5}));
    auto x_plus_1 = cc->EvalAdd(x, homomorphic_1);
    auto result = cc->EvalMult(x_plus_1, homomorphic_05);
    return result;
}


vector<vector<Ciphertext<DCRTPoly>>> compute_comparisons(const CryptoContext<DCRTPoly> &cc, const vector<vector<Ciphertext<DCRTPoly>>> &A, const Ciphertext<DCRTPoly> &c_enc_05){

    //Number of compositions
    int df=2;
    int dg=2;
    int vec_size=1;

    vector<vector<Ciphertext<DCRTPoly>>> comparisonSetA(A.size(), vector<Ciphertext<DCRTPoly>>(A.size()));
    for (size_t i=0; i<A.size(); i++){
        for (size_t j=i; j<A.size(); j++){
            if (i!=j){
                comparisonSetA[i][j] = homomorphic_comparison_g(cc,A[i][0],A[j][0],df,dg,vec_size);
                comparisonSetA[j][i] = cc->EvalAdd(1,-comparisonSetA[i][j]);
            } else {
                comparisonSetA[i][j] = c_enc_05;
            }
        }
    }
    return comparisonSetA;
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
     vector<vector<double>> dataset = {
        {50.2, 3.6, 1.9, 1},
        {21.8, 1.7, 0.8, 0},
        {40.5, 2.3, 0.9, 1},
        {35.3, 3.2, 1.4, 0},
        {32.0, 2.2, 1.3, 1},
        {45.3, 3.1, 1.7, 1},
        {38.9, 2.8, 1.4, 1},
        {29.1, 1.9, 0.9, 0},
        {25.5, 1.6, 0.8, 0},

        {50.2, 3.6, 1.9, 1},
        {21.8, 1.7, 0.8, 0},
        {40.5, 2.3, 0.9, 1},
        {35.3, 3.2, 1.4, 0},
        {32.0, 2.2, 1.3, 1},
        {45.3, 3.1, 1.7, 1},
        {38.9, 2.8, 1.4, 1},
        {29.1, 1.9, 0.9, 0},
        {25.5, 1.6, 0.8, 0},
    };

    vector<vector<vector<double>>> new_dataset;
    transform(dataset.begin(), dataset.end(), back_inserter(new_dataset), 
              [](const vector<double>& row) {
                  vector<vector<double>> row_vectors;
                  transform(row.begin(), row.end(), back_inserter(row_vectors), 
                            [](double x) { return vector<double>{x}; });
                  return row_vectors;
              });

    for (const auto& row : new_dataset) {
        for (const auto& element : row) {
            cout << "{" << element[0] << "}, ";
        }
        cout << endl;
    }
    cout << endl;

    
    //STEP 1 - Escalar dataset
    auto scaled_dataset = min_max_scaling(new_dataset);

    for (const auto& row : scaled_dataset) {
        for (const auto& element : row) {
            cout << "{" << element[0] << "}, ";
        }
        cout << endl;
    }
    cout << endl;


    //STEP 2- Encriptar dades
    auto global_start = std::chrono::high_resolution_clock::now();


    //2.1 Setup Parameters FHE
    uint32_t multDepth = 40;
    uint32_t scaleModSize = 50;
    uint32_t batchSize = 1;

    CCParams<CryptoContextCKKSRNS> parameters;

    parameters.SetMultiplicativeDepth(multDepth);
    parameters.SetScalingModSize(scaleModSize);
    parameters.SetScalingTechnique(FLEXIBLEAUTO);
    parameters.SetBatchSize(batchSize);

    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);

    auto keys = cc->KeyGen();
    cc->EvalMultKeyGen(keys.secretKey);

    //2.2 Encrypt
    vector<vector<Ciphertext<DCRTPoly>>> encrypted_dataset;
    for (const auto& row : scaled_dataset) {
        vector<Ciphertext<DCRTPoly>> encrypted_row(row.size());
        for (size_t i = 0; i < row.size(); i++) {
            vector<double> values = { row[i][0] };
            Plaintext ptxt = cc->MakeCKKSPackedPlaintext(values);
            encrypted_row[i] = cc->Encrypt(keys.publicKey, ptxt);
        }
        encrypted_dataset.push_back(encrypted_row);
    }

    //2.EXTRA - Desencriptar dades
    /*vector<vector<Plaintext>> decrypted_dataset;
    for (const auto& row : encrypted_dataset) {
        vector<Plaintext> decrypted_row(row.size());
        for (size_t i = 0; i < row.size(); i++) {
            Plaintext plaintext;
            cc->Decrypt(keys.secretKey, row[i], &plaintext);
            decrypted_row[i] = plaintext;
        }
        decrypted_dataset.push_back(decrypted_row);
    }

    // Print the decrypted dataset
    for (const auto& row : decrypted_dataset) {
        for (const auto& element : row) {
            cout << "{" << element << "}, ";
        }
        cout << endl;
    }*/

    //2.3 Separar train i test
    
    //2.3.1 De forma aleatoria
    /*vector<vector<Ciphertext<DCRTPoly>>> shuffled_dataset = encrypted_dataset;
    random_device rd;
    mt19937 g(rd());
    shuffle(shuffled_dataset.begin(), shuffled_dataset.end(), g);

    int num_samples = shuffled_dataset.size();
    int num_train_samples = 0.9 * num_samples;
    vector<vector<Ciphertext<DCRTPoly>>> train_dataset(shuffled_dataset.begin(), shuffled_dataset.begin() + num_train_samples);
    vector<vector<Ciphertext<DCRTPoly>>> test_dataset(shuffled_dataset.begin() + num_train_samples, shuffled_dataset.end());*/

    //2.3.2 Els primers son train i els ultims son test
    int num_samples = encrypted_dataset.size();
    int num_train_samples = 0.99 * num_samples;
    vector<vector<Ciphertext<DCRTPoly>>> train_dataset(encrypted_dataset.begin(), encrypted_dataset.begin() + num_train_samples);
    vector<vector<Ciphertext<DCRTPoly>>> test_dataset(encrypted_dataset.begin() + num_train_samples, encrypted_dataset.end());

    cout << "Train samples: " << train_dataset.size() << endl;
    cout << "Test samples: " << test_dataset.size() << endl;
    cout << endl;

    //Encryption of 0.5. Variable auxiliar que es fa servir al pas 4.1
    std::vector<double> enc_05 = {0.5};
    Plaintext p_05 = cc->MakeCKKSPackedPlaintext(enc_05);
    auto c_enc_05 = cc->Encrypt(keys.publicKey, p_05);


    //Variable auxiliar que es fa servir al pas 4.2
    std::vector<double> aux = {-1.0};
    Plaintext p_aux = cc->MakeCKKSPackedPlaintext(aux);
    auto c_aux = cc->Encrypt(keys.publicKey, p_aux);


    for (size_t test_index=0; test_index<test_dataset.size(); test_index++){
        
        //Agafar els features, no la label
        auto row_test = test_dataset[test_index];
        vector<Ciphertext<DCRTPoly>> features_test(row_test.begin(), row_test.end() - 1);

        //Variable (array de distancies i labels) que s'ha d'ordenar
        vector<vector<Ciphertext<DCRTPoly>>> A;

        for (size_t train_index=0; train_index<train_dataset.size(); train_index++){
            auto row_train = train_dataset[train_index];
            vector<Ciphertext<DCRTPoly>> features_train(row_train.begin(), row_train.end() - 1);
            Ciphertext<DCRTPoly> label_train = row_train[row_train.size()-1];
            //STEP 3 - Calcular euclidean distance (per cada sample de test)
            Ciphertext<DCRTPoly> dist = euclidean_distance(cc,features_train,features_test);

            vector<Ciphertext<DCRTPoly>> dist_and_label = {dist, label_train};
            A.push_back(dist_and_label);
        }

        //STEP 4 - Ordenar euclidean distance (i respectiva label), fent servir merge_sort
        //4.1 Calcular pairwise comparisons
        cout << "Computing comparisons..." << endl;

        auto start = std::chrono::high_resolution_clock::now();
        vector<vector<Ciphertext<DCRTPoly>>> compSetA = compute_comparisons(cc,A,c_enc_05);
        auto stop_comp = std::chrono::high_resolution_clock::now();
        auto duration_comp = std::chrono::duration_cast<std::chrono::microseconds>(stop_comp - start);

        cout << "Comparisons computed" << endl;
        std::cout << "Time taken in computing comparisons: " << std::chrono::duration<double>(duration_comp).count() << " seconds" << std::endl;
        cout << endl;

        //4.2 Ordenar (Merge_sort)
        cout << "Ordering values..." << endl;

        vector<vector<Ciphertext<DCRTPoly>>> sorted_A = merge_sort(cc,A,compSetA,c_aux);
        auto stop_merge = std::chrono::high_resolution_clock::now();
        auto duration_merge = std::chrono::duration_cast<std::chrono::microseconds>(stop_merge - stop_comp);

        cout << "Distances ordered" << endl;
        std::cout << "Time taken merge_sort: " << std::chrono::duration<double>(duration_merge).count() << " seconds" << std::endl;
        cout << endl;



        //STEP 5 - Agafar labels de les k distàncies més petites
        int k=2;
        vector<vector<Ciphertext<DCRTPoly>>> last_k_elements(sorted_A.end()-k,sorted_A.end());
        vector<Ciphertext<DCRTPoly>> last_k_labels;
        for (size_t i=0; i<last_k_elements.size(); i++){
            auto label_i = last_k_elements[i][1];
            last_k_labels.push_back(label_i);
        }


        //STEP 6 - Sumar les k labels
        auto sum_labels = last_k_labels[0];
        for (size_t i=1; i<last_k_labels.size(); i++){
            sum_labels = cc->EvalAdd(sum_labels,last_k_labels[i]);
        }


        //STEP 7 - Desencriptar resultat de step 6) i en funció d'això decidir la label
        Plaintext result;
        cc->Decrypt(keys.secretKey, sum_labels, &result);
        cout << "Sum of labels: " << result << endl;
        
    }

    auto global_stop = std::chrono::high_resolution_clock::now();
    auto total_duration = std::chrono::duration_cast<std::chrono::microseconds>(global_stop - global_start);
    std::cout << "Total time taken: " << std::chrono::duration<double>(total_duration).count() << " seconds" << std::endl;

    return 0;
}