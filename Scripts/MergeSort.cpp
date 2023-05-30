#define PROFILE

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "openfhe.h"
#include <cmath>

using namespace lbcrypto;
using namespace std;


Ciphertext<DCRTPoly> f_4_homomorphic(const CryptoContext<DCRTPoly> &cc, const Ciphertext<DCRTPoly> &c) {
    auto c2 = cc->EvalMult(c, c); // c^2
    auto c3 = cc->EvalMult(c2, c); // c^3
    auto c5 = cc->EvalMult(c2, c3); // c^5
    auto c7 = cc->EvalMult(c2, c5); // c^7
    auto c9 = cc->EvalMult(c2, c7); // c^9

    auto term1 = cc->EvalMult(35.0/128.0, c9);
    auto term2 = cc->EvalMult(-180.0/128.0, c7);
    auto term3 = cc->EvalMult(378.0/128.0, c5);
    auto term4 = cc->EvalMult(-420.0/128.0, c3);
    auto term5 = cc->EvalMult(315.0/128.0, c);

    auto result = cc->EvalAdd(term1, term2);
    result = cc->EvalAdd(result, term3);
    result = cc->EvalAdd(result, term4);
    result = cc->EvalAdd(result, term5);

    return result;
}

Ciphertext<DCRTPoly> g_4_homomorphic(const CryptoContext<DCRTPoly> &cc, const Ciphertext<DCRTPoly> &c) {
    auto c2 = cc->EvalMult(c, c); // c^2
    auto c3 = cc->EvalMult(c2, c); // c^3
    auto c5 = cc->EvalMult(c2, c3); // c^5
    auto c7 = cc->EvalMult(c2, c5); // c^7
    auto c9 = cc->EvalMult(c2, c7); // c^9

    auto term1 = cc->EvalMult(46623.0 / 1024.0, c9);
    auto term2 = cc->EvalMult(-113492.0 / 1024.0, c7);
    auto term3 = cc->EvalMult(97015.0 / 1024.0, c5);
    auto term4 = cc->EvalMult(-34974.0 / 1024.0, c3);
    auto term5 = cc->EvalMult(5850.0 / 1024.0, c);

    auto result = cc->EvalAdd(term1, term2);
    result = cc->EvalAdd(result, term3);
    result = cc->EvalAdd(result, term4);
    result = cc->EvalAdd(result, term5);

    return result;
}


Ciphertext<DCRTPoly> homomorphic_comparison_g(const CryptoContext<DCRTPoly> &cc, const Ciphertext<DCRTPoly> &c1, const Ciphertext<DCRTPoly> &c2) {
    int df=2;
    int dg=2;
    
    auto x = cc->EvalSub(c1, c2);
    for (int i=1; i<=dg; i++){
        x = g_4_homomorphic(cc,x);
    }
    for (int i=1; i<=df; i++){
        x = f_4_homomorphic(cc,x);
    }

    auto x_plus_1 = cc->EvalAdd(x, 1);
    auto result = cc->EvalMult(x_plus_1, 0.5);
    return result;
}


vector<Ciphertext<DCRTPoly>> L_function(const CryptoContext<DCRTPoly> &cc, const Ciphertext<DCRTPoly> comparison, const vector<Ciphertext<DCRTPoly>> &F, const vector<Ciphertext<DCRTPoly>> &G ){

    auto invComp = cc->EvalAdd(1,-comparison);

    //L function calculation for the first element of the tuple (distance)
    auto first_term_dist = cc->EvalMult(comparison,F[0]);
    auto second_term_dist = cc->EvalMult(invComp,G[0]);
    auto result_dist = cc->EvalAdd(first_term_dist,second_term_dist);

    //L function calculation for the second element of the tuple (label)
    auto first_term_label = cc->EvalMult(comparison,F[1]);
    auto second_term_label = cc->EvalMult(invComp,G[1]);
    auto result_label = cc->EvalAdd(first_term_label,second_term_label);

    vector<Ciphertext<DCRTPoly>> result = {result_dist,result_label};
    return result;
}

vector<Ciphertext<DCRTPoly>> m_max(int m, const CryptoContext<DCRTPoly> &cc, const vector<vector<Ciphertext<DCRTPoly>>> &B, const vector<vector<Ciphertext<DCRTPoly>>> &C,const vector<vector<Ciphertext<DCRTPoly>>> &comparisonSetBC) 
{   
    int size_b = B.size();
    int size_c = C.size();

    //base case 1: one array is empty
    if (size_b == 0 || size_c == 0) {
       return size_b == 0 ? C[m - 1] : B[m - 1];
    } 

    //base case 2: when M=1 return max of first elements of both arrays
    if (m==1) {
        return L_function(cc, comparisonSetBC[0][0], B[0], C[0]);
    } 
    
    //base case 3: when M=size_b+size_c return min of last elements of both arrays
    if (m==size_b+size_c){
        return L_function(cc, comparisonSetBC[size_b-1][size_c-1], C[size_c-1], B[size_b-1]);
    }


    int i = floor(m / 2.0);
    int j = ceil(m / 2.0);

    vector<Ciphertext<DCRTPoly>> left;
    vector<Ciphertext<DCRTPoly>> right;

    int n=comparisonSetBC.size();
    int p=comparisonSetBC[0].size();

    if (i<size_b){ 
        if (j <= size_c) {
            #pragma omp task shared(left)
            {
            vector<vector<Ciphertext<DCRTPoly>>> comparisonSetLEFT(n - i);
            #pragma omp parallel for
                for (int x = 0; x < n - i; x++) {
                    comparisonSetLEFT[x].resize(j);
                    for (int y = 0; y < j; y++) {
                        comparisonSetLEFT[x][y] = comparisonSetBC[i + x][y];
                    }
                }
                left = m_max(j, cc, vector<vector<Ciphertext<DCRTPoly>>>(B.begin() + i, B.end()), vector<vector<Ciphertext<DCRTPoly>>>(C.begin(), C.begin() + j), comparisonSetLEFT);
            }

            #pragma omp task shared(right)
            {
            vector<vector<Ciphertext<DCRTPoly>>> comparisonSetRIGHT(i);
            #pragma omp parallel for
                for (int x = 0; x < i; x++) {
                    comparisonSetRIGHT[x].resize(p - j);
                    for (int y = 0; y < p - j; y++) {
                        comparisonSetRIGHT[x][y] = comparisonSetBC[x][j + y];
                    }
                }
                right = m_max(i, cc, vector<vector<Ciphertext<DCRTPoly>>>(B.begin(), B.begin() + i), vector<vector<Ciphertext<DCRTPoly>>>(C.begin() + j, C.end()), comparisonSetRIGHT);
            }
            } 
        
        else 
        {
            vector<vector<Ciphertext<DCRTPoly>>> comparisonSetLEFT(n-i);
            #pragma omp parallel for
            for (int x = 0; x < n-i; x++) {
                comparisonSetLEFT[x].resize(p);
                for (int y = 0; y < p; y++) {
                    comparisonSetLEFT[x][y] = comparisonSetBC[i+x][y];
                }
            }
            #pragma omp task shared(left)
            {
            left = m_max(j, cc, vector<vector<Ciphertext<DCRTPoly>>>(B.begin() + i, B.end()), vector<vector<Ciphertext<DCRTPoly>>>(C.begin(), C.end()),comparisonSetLEFT);
            }


            vector<vector<Ciphertext<DCRTPoly>>> comparisonSetRIGHT = vector<vector<Ciphertext<DCRTPoly>>>(comparisonSetBC.begin(), comparisonSetBC.begin());
            #pragma omp task shared(right)
            {
            right = m_max(i, cc, vector<vector<Ciphertext<DCRTPoly>>>(B.begin(), B.begin() + i), vector<vector<Ciphertext<DCRTPoly>>>(C.begin(), C.begin()),comparisonSetRIGHT);
            }
        } 
    }
    else { // i>size_b

        vector<vector<Ciphertext<DCRTPoly>>> comparisonSetLEFT = vector<vector<Ciphertext<DCRTPoly>>>(comparisonSetBC.begin(), comparisonSetBC.begin());
        #pragma omp task shared(left)
        {
        left = m_max(j, cc, vector<vector<Ciphertext<DCRTPoly>>>(B.begin(), B.begin()), vector<vector<Ciphertext<DCRTPoly>>>(C.begin(), C.begin() + j),comparisonSetLEFT);
        }

        vector<vector<Ciphertext<DCRTPoly>>> comparisonSetRIGHT(n);
            #pragma omp parallel for
            for (int x = 0; x < n; x++) {
                comparisonSetRIGHT[x].resize(p - j);
                for (int y = 0; y < p-j; y++) {
                    comparisonSetRIGHT[x][y] = comparisonSetBC[x][j+y];
                }
            }
        #pragma omp task shared(right)
        {
        right = m_max(i, cc, vector<vector<Ciphertext<DCRTPoly>>>(B.begin(), B.end()), vector<vector<Ciphertext<DCRTPoly>>>(C.begin() + j, C.end()),comparisonSetRIGHT);
        }
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
    #pragma omp taskwait
    return L_function(cc, xiCOMPyj, left, right);
}



vector<vector<Ciphertext<DCRTPoly>>> merge(const CryptoContext<DCRTPoly> &cc, const Ciphertext<DCRTPoly> c0, const vector<vector<Ciphertext<DCRTPoly>>> &B, const vector<vector<Ciphertext<DCRTPoly>>> &C, const vector<vector<Ciphertext<DCRTPoly>>> &comparisonSetBC){
    int s = B.size();
    int t = C.size();

    
    vector<vector<Ciphertext<DCRTPoly>>> Z(s+t);

    //1) s+t-1 calls to m_max function
    #pragma omp parallel for schedule(dynamic)
    for (int i=0; i<(s+t-2); i++){
        vector<Ciphertext<DCRTPoly>> zi = m_max(i+1,cc,B,C,comparisonSetBC);
        Z[i] = zi;
        
    }
    Z[s+t-1] = m_max(s+t,cc,B,C,comparisonSetBC);

    //2) Compute missing element

    //sum all elements B
    Ciphertext<DCRTPoly> zi_dist=c0;
    Ciphertext<DCRTPoly> zi_label=c0;
    for (int i=0; i<s; i++){
        zi_dist = cc->EvalAdd(zi_dist,B[i][0]);
        zi_label = cc->EvalAdd(zi_label,B[i][1]);
    }
    // sum all elements in C
    for (int i=0; i<t; i++){
        zi_dist = cc->EvalAdd(zi_dist,C[i][0]);
        zi_label = cc->EvalAdd(zi_label,C[i][1]);
    }

    //Subtract all elements computed in Z (except last element)
    for (int i=0; i<(s+t-2); i++){
        zi_dist = cc->EvalSub(zi_dist,Z[i][0]);
        zi_label = cc->EvalSub(zi_label,Z[i][1]);
    }
    zi_dist = cc->EvalSub(zi_dist, Z[s+t-1][0]);
    zi_label = cc->EvalSub(zi_label, Z[s+t-1][1]);

    Z[s+t-2]={zi_dist,zi_label};

    return Z;
}

vector<vector<Ciphertext<DCRTPoly>>> merge_sort(const CryptoContext<DCRTPoly> &cc, const Ciphertext<DCRTPoly> &c0, const vector<vector<Ciphertext<DCRTPoly>>> &A, const vector<vector<Ciphertext<DCRTPoly>>> &comparisonSetA){
    
    //Base case
    if (A.size()==1) {
        return A;
    }
    int s=floor(A.size()/2);
    int n=A.size(); 
    
    //Step 1: Recursively sort left and right parts of the array
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

    vector<vector<Ciphertext<DCRTPoly>>> B, C;

    #pragma omp task shared(B)
        {  
            B = merge_sort(cc,c0,vector<vector<Ciphertext<DCRTPoly>>>(A.begin(), A.begin()+s), comparisonSetB); //Recursively sort left part of the array
        }

    #pragma omp task shared(C)
        {
            C = merge_sort(cc,c0,vector<vector<Ciphertext<DCRTPoly>>>(A.begin()+s, A.end()), comparisonSetC); //Recursively sort right part of the array
        } 

    int k=A.size();

    //Step 2: Sort comparison results
    vector<vector<vector<Ciphertext<DCRTPoly>>>> unionBiAj(k-s);
    auto sorterParameter = vector<vector<Ciphertext<DCRTPoly>>>(comparisonSetA.begin(), comparisonSetA.begin()+s);
    #pragma omp parallel for
    for (int j = (s+1); j<=k; j++){
    
        vector<vector<Ciphertext<DCRTPoly>>> sorterParameter1;

        for (int i=0; i<s; i++){
            vector<Ciphertext<DCRTPoly>> row = {comparisonSetA[i][j-1], c0};
            sorterParameter1.push_back(row);
        }

        auto BCOMPaj = merge_sort(cc,c0,sorterParameter1,comparisonSetB);
        unionBiAj[j-s-1] = BCOMPaj;
    }



    vector<vector<vector<Ciphertext<DCRTPoly>>>> unionBC(s);
    #pragma omp parallel for
    for (int i = 1; i<=s; i++){
        vector<vector<Ciphertext<DCRTPoly>>> sorterParameter2;
        for(int j=0; j<k-s; j++){ 
            sorterParameter2.push_back(unionBiAj[j][i-1]);
        }
        auto biCOMPC = merge_sort(cc,c0,sorterParameter2,comparisonSetC);
        unionBC[i-1] = biCOMPC;
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
    #pragma omp taskwait

    //Step 3: Merge the two sorted halves using the new comparison matrix
    auto return_value = merge(cc,c0,B,C,BCompC);
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


    vector<vector<Ciphertext<DCRTPoly>>> A = {c5,c3,c2,c6,c1,c8,c7,c4};

    //Compute pairwise comparisons
    std::vector<double> enc_05 = {0.5};
    Plaintext p_05 = cc->MakeCKKSPackedPlaintext(enc_05);
    auto c_enc_05 = cc->Encrypt(keys.publicKey, p_05);

    auto start = std::chrono::high_resolution_clock::now();


    vector<vector<Ciphertext<DCRTPoly>>> comparisonSetA(8, vector<Ciphertext<DCRTPoly>>(8));
    for (size_t i=0; i<A.size(); i++){
        for (size_t j=i; j<A.size(); j++){
            if (i!=j){
            comparisonSetA[i][j] = homomorphic_comparison_g(cc,A[i][0],A[j][0]);
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


    std::vector<double> aux = {0.0};
    Plaintext p_aux = cc->MakeCKKSPackedPlaintext(aux);
    auto c_aux = cc->Encrypt(keys.publicKey, p_aux);
    vector<vector<Ciphertext<DCRTPoly>>> sorted_A = merge_sort(cc,c_aux,A,comparisonSetA);

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
    
    


    

    
   

    
