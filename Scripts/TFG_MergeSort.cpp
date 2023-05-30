#define PROFILE

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "openfhe.h"
#include <cmath>

using namespace lbcrypto;
using namespace std;


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
        auto sorterParameter = vector<vector<Ciphertext<DCRTPoly>>>(comparisonSetA.begin(), comparisonSetA.begin()+s); 
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
    
    // distance, label
    std::vector<double> distance_1 = {0.1};
    std::vector<double> label_1 = {1.0};

    std::vector<double> distance_2 = {0.2};
    std::vector<double> label_2 = {0.0};

    std::vector<double> distance_3 = {0.3};
    std::vector<double> label_3 = {1.0};

    std::vector<double> distance_4 = {0.4};
    std::vector<double> label_4 = {0.0};

    std::vector<double> distance_5 = {0.5};
    std::vector<double> label_5 = {0.0};

    std::vector<double> distance_6 = {0.6};
    std::vector<double> label_6 = {0.0};

    std::vector<double> distance_7 = {0.7};
    std::vector<double> label_7 = {1.0};

    std::vector<double> distance_8 = {0.8};
    std::vector<double> label_8 = {0.0};

    // Setup parameters cryptocontext
    uint32_t multDepth = 38;
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


    Plaintext ptxt1 = cc->MakeCKKSPackedPlaintext(distance_1);
    Plaintext ptxt2 = cc->MakeCKKSPackedPlaintext(distance_2);
    Plaintext ptxt3 = cc->MakeCKKSPackedPlaintext(distance_3);
    Plaintext ptxt4 = cc->MakeCKKSPackedPlaintext(distance_4);
    Plaintext ptxt5 = cc->MakeCKKSPackedPlaintext(distance_5);
    Plaintext ptxt6 = cc->MakeCKKSPackedPlaintext(distance_6);
    Plaintext ptxt7 = cc->MakeCKKSPackedPlaintext(distance_7);
    Plaintext ptxt8 = cc->MakeCKKSPackedPlaintext(distance_8);

    Plaintext ptxtlabel1 = cc->MakeCKKSPackedPlaintext(label_1);
    Plaintext ptxtlabel2 = cc->MakeCKKSPackedPlaintext(label_2);
    Plaintext ptxtlabel3 = cc->MakeCKKSPackedPlaintext(label_3);
    Plaintext ptxtlabel4 = cc->MakeCKKSPackedPlaintext(label_4);
    Plaintext ptxtlabel5 = cc->MakeCKKSPackedPlaintext(label_5);
    Plaintext ptxtlabel6 = cc->MakeCKKSPackedPlaintext(label_6);
    Plaintext ptxtlabel7 = cc->MakeCKKSPackedPlaintext(label_7);
    Plaintext ptxtlabel8 = cc->MakeCKKSPackedPlaintext(label_8);


    // Encrypt
    auto cd1 = cc->Encrypt(keys.publicKey, ptxt1);
    auto cd2 = cc->Encrypt(keys.publicKey, ptxt2);
    auto cd3 = cc->Encrypt(keys.publicKey, ptxt3);
    auto cd4 = cc->Encrypt(keys.publicKey, ptxt4);
    auto cd5 = cc->Encrypt(keys.publicKey, ptxt5);
    auto cd6 = cc->Encrypt(keys.publicKey, ptxt6);
    auto cd7 = cc->Encrypt(keys.publicKey, ptxt7);
    auto cd8 = cc->Encrypt(keys.publicKey, ptxt8);

    auto cl1 = cc->Encrypt(keys.publicKey, ptxtlabel1);
    auto cl2 = cc->Encrypt(keys.publicKey, ptxtlabel2);
    auto cl3 = cc->Encrypt(keys.publicKey, ptxtlabel3);
    auto cl4 = cc->Encrypt(keys.publicKey, ptxtlabel4);
    auto cl5 = cc->Encrypt(keys.publicKey, ptxtlabel5);
    auto cl6 = cc->Encrypt(keys.publicKey, ptxtlabel6);
    auto cl7 = cc->Encrypt(keys.publicKey, ptxtlabel7);
    auto cl8 = cc->Encrypt(keys.publicKey, ptxtlabel8);

    vector<Ciphertext<DCRTPoly>> c1 = {cd1, cl1};
    vector<Ciphertext<DCRTPoly>> c2 = {cd2, cl2};
    vector<Ciphertext<DCRTPoly>> c3 = {cd3, cl3};
    vector<Ciphertext<DCRTPoly>> c4 = {cd4, cl4};
    vector<Ciphertext<DCRTPoly>> c5 = {cd5, cl5};
    vector<Ciphertext<DCRTPoly>> c6 = {cd6, cl6};
    vector<Ciphertext<DCRTPoly>> c7 = {cd7, cl7};
    vector<Ciphertext<DCRTPoly>> c8 = {cd8, cl8};


    vector<vector<Ciphertext<DCRTPoly>>> A = {c5,c3,c2,c6,c1,c8,c5,c4};

    int df = 2;
    int dg = 2;
    int vec_size = 1;

    //Compute pairwise comparisons
    std::vector<double> enc_05 = {0.5};
    Plaintext p_05 = cc->MakeCKKSPackedPlaintext(enc_05);
    auto c_enc_05 = cc->Encrypt(keys.publicKey, p_05);



    auto start = std::chrono::high_resolution_clock::now();


    vector<vector<Ciphertext<DCRTPoly>>> comparisonSetA(8, vector<Ciphertext<DCRTPoly>>(8));
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

    cout << "comparisons computed" << endl;


    auto stop_comp = std::chrono::high_resolution_clock::now();
    auto duration_comp = std::chrono::duration_cast<std::chrono::microseconds>(stop_comp - start);
    std::cout << "Time taken in computing comparisons: " << std::chrono::duration<double>(duration_comp).count() << " seconds" << std::endl;


    std::vector<double> aux = {-1.0};
    Plaintext p_aux = cc->MakeCKKSPackedPlaintext(aux);
    auto c_aux = cc->Encrypt(keys.publicKey, p_aux);
    vector<vector<Ciphertext<DCRTPoly>>> sorted_A = merge_sort(cc,A,comparisonSetA,c_aux);

    auto stop_merge = std::chrono::high_resolution_clock::now();
    auto duration_merge = std::chrono::duration_cast<std::chrono::microseconds>(stop_merge - stop_comp);
    std::cout << "Time taken merge_sort: " << std::chrono::duration<double>(duration_merge).count() << " seconds" << std::endl;


    auto total_duration = std::chrono::duration_cast<std::chrono::microseconds>(stop_merge - start);
    std::cout << "Total time taken: " << std::chrono::duration<double>(total_duration).count() << " seconds" << std::endl;


    for (size_t i=0; i<A.size(); i++){
        Plaintext plaintext_dist;
        Plaintext plaintext_label;
        cc->Decrypt(keys.secretKey, sorted_A[i][0], &plaintext_dist);
        cc->Decrypt(keys.secretKey, sorted_A[i][1], &plaintext_label);
        cout << "Decrypted value (distance) of " << i << "-th max element: " << plaintext_dist << endl;
        cout << "Decrypted value (label) of " << i << "-th max element: " << plaintext_label << endl;
        }

    
    return 0;

}
    
    


    

    
   

    
