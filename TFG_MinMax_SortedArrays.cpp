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

/*vector<Ciphertext<DCRTPoly>> homomorphic_max(const CryptoContext<DCRTPoly> &cc, const vector<Ciphertext<DCRTPoly>> &ca, const vector<Ciphertext<DCRTPoly>> &cb, const Ciphertext<DCRTPoly> &acompb)
{

    auto first_term_dist = cc->EvalMult(acompb,ca[0]);
    auto second_term_dist = cc->EvalMult(cc->EvalAdd(1,-acompb),cb[0]);
    auto result_dist = cc->EvalAdd(first_term_dist,second_term_dist);

    auto first_term_label = cc->EvalMult(acompb,ca[1]);
    auto second_term_label = cc->EvalMult(cc->EvalAdd(1,-acompb),cb[1]);
    auto result_label = cc->EvalAdd(first_term_label,second_term_label);

    vector<Ciphertext<DCRTPoly>> result = {result_dist,result_label};
    return result;
}*/

/*vector<Ciphertext<DCRTPoly>> homomorphic_min(const CryptoContext<DCRTPoly> &cc, const vector<Ciphertext<DCRTPoly>> &ca, const vector<Ciphertext<DCRTPoly>> &cb, const Ciphertext<DCRTPoly> &acompb)
// min(a, b) = (a < b) * a + (a > b) * b
{

    auto first_term_dist = cc->EvalMult(cc->EvalAdd(1,-acompb),ca[0]);
    auto second_term_dist = cc->EvalMult(acompb,cb[0]);
    auto result_dist = cc->EvalAdd(first_term_dist,second_term_dist);

    auto first_term_label = cc->EvalMult(cc->EvalAdd(1,-acompb),ca[1]);
    auto second_term_label = cc->EvalMult(acompb,cb[1]);
    auto result_label = cc->EvalAdd(first_term_label,second_term_label);

    vector<Ciphertext<DCRTPoly>> result = {result_dist,result_label};
    return result;
}*/


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



vector<Ciphertext<DCRTPoly>> m_max(int m, const CryptoContext<DCRTPoly> &cc, const vector<vector<Ciphertext<DCRTPoly>>> &B, const vector<vector<Ciphertext<DCRTPoly>>> &C,const vector<vector<Ciphertext<DCRTPoly>>>> &comparisonSetBC) 
{   
    int df = 2;
    int dg = 2;
    int vec_size = 1;

    int size_b = B.size();
    int size_c = C.size();

    //base case
    if (size_b == 0 || size_c == 0) {
       return size_b == 0 ? C[m - 1] : B[m - 1];
    } 

    //base case, when M=1 return max of first elements of both arrays
    if (m==1) {
        auto xiCOMPyj = homomorphic_comparison_g(cc, B[0][0], C[0][0], df, dg, vec_size);
        return L_function(cc, xiCOMPyj, B[0], C[0]);
    } 
    
    //base case, when M=size_b+size_c return min of last elements of both arrays
    if (m==size_b+size_c){
        auto xiCOMPyj = homomorphic_comparison_g(cc, B[size_b-1][0], C[size_c-1][0], df, dg, vec_size);
        return L_function(cc, xiCOMPyj, C[size_c-1], B[size_b-1]);
    }


    int i = floor(m / 2.0);
    int j = ceil(m / 2.0);

    vector<Ciphertext<DCRTPoly>> left;
    vector<Ciphertext<DCRTPoly>> right;

    if (i<=size_b){
        if (j<=size_c){
            left = m_max(j, cc, vector<vector<Ciphertext<DCRTPoly>>>(B.begin() + i, B.end()), vector<vector<Ciphertext<DCRTPoly>>>(C.begin(), C.begin() + j));
            right = m_max(i, cc, vector<vector<Ciphertext<DCRTPoly>>>(B.begin(), B.begin() + i), vector<vector<Ciphertext<DCRTPoly>>>(C.begin() + j, C.end()));
        } else {
            left = m_max(j, cc, vector<vector<Ciphertext<DCRTPoly>>>(B.begin() + i, B.end()), vector<vector<Ciphertext<DCRTPoly>>>(C.begin(), C.end()));
            right = m_max(i, cc, vector<vector<Ciphertext<DCRTPoly>>>(B.begin(), B.begin() + i), vector<vector<Ciphertext<DCRTPoly>>>(C.begin(), C.begin()));
        } 
    }
    else {
        left = m_max(j, cc, vector<vector<Ciphertext<DCRTPoly>>>(B.begin(), B.begin()), vector<vector<Ciphertext<DCRTPoly>>>(C.begin(), C.begin() + j));
        right = m_max(i, cc, vector<vector<Ciphertext<DCRTPoly>>>(B.begin(), B.end()), vector<vector<Ciphertext<DCRTPoly>>>(C.begin() + j, C.end()));
    }

    Ciphertext<DCRTPoly> xiCOMPyj;
    if (i>size_b){
        xiCOMPyj = homomorphic_comparison_g(cc, B[size_b-1][0], C[j-1][0], df, dg, vec_size);
    } else {
        if (j>size_c){
            xiCOMPyj = homomorphic_comparison_g(cc, B[i-1][0], C[size_c-1][0],df,dg,vec_size);
        }
        else{
            xiCOMPyj = homomorphic_comparison_g(cc, B[i-1][0], C[j-1][0], df, dg, vec_size);
        }
    }
    return L_function(cc, xiCOMPyj, left, right);
}




/*vector<Ciphertext<DCRTPoly>> m_min(int m, const CryptoContext<DCRTPoly> &cc, const vector<vector<Ciphertext<DCRTPoly>>> &B, const vector<vector<Ciphertext<DCRTPoly>>> &C) //const vector<vector<Ciphertext<DCRTPoly>>> &BCompC
{   
    int df = 2;
    int dg = 2;
    int vec_size = 1;

    int size_b = B.size();
    int size_c = C.size();
    int s = size_b ;
    int t = size_c ;

    if (size_b == 0 || size_c == 0) {
        if (size_b == 0) {
            return C[t-m];
        } else {
            return B[s-m];
        }
    } else {
        int i = s-floor(m / 2.0);
        int j = t-ceil(m / 2.0);


        vector<Ciphertext<DCRTPoly>> left;
        if (i>=1){
            left = m_min(j, cc, vector<vector<Ciphertext<DCRTPoly>>>(B.begin() , B.begin()+i-1), vector<vector<Ciphertext<DCRTPoly>>>(C.begin()+j, C.end()));
        } else{
             left = m_min(j, cc, vector<vector<Ciphertext<DCRTPoly>>>(B.begin() , B.begin()), vector<vector<Ciphertext<DCRTPoly>>>(C.begin()+j, C.end()));
        }

        vector<Ciphertext<DCRTPoly>> right;
        if (j>=1){
            right = m_min(i, cc, vector<vector<Ciphertext<DCRTPoly>>>(B.begin()+i, B.end()), vector<vector<Ciphertext<DCRTPoly>>>(C.begin(), C.begin()+j-1));    
        } else {
            right = m_min(i, cc, vector<vector<Ciphertext<DCRTPoly>>>(B.begin()+i, B.end()), vector<vector<Ciphertext<DCRTPoly>>>(C.begin(), C.begin()));
        }
        auto leftCOMPright = homomorphic_comparison_g(cc, left[0], right[0], df, dg, vec_size);

        return homomorphic_min(cc, left, right, leftCOMPright);
    }
}*/



int main() {
    
    /*
    int df=2;
    int dg=2;
    int vec_size=1;*/

    // distance, label
    std::vector<double> distance_1 = {0.2};
    std::vector<double> label_1 = {1.0};

    std::vector<double> distance_2 = {0.3};
    std::vector<double> label_2 = {0.0};

    std::vector<double> distance_3 = {0.4};
    std::vector<double> label_3 = {1.0};

    std::vector<double> distance_4 = {0.5};
    std::vector<double> label_4 = {0.0};

    std::vector<double> distance_5 = {0.4};
    std::vector<double> label_5 = {0.0};

    std::vector<double> distance_6 = {0.7};
    std::vector<double> label_6 = {0.0};

    std::vector<double> distance_7 = {0.9};
    std::vector<double> label_7 = {1.0};

    std::vector<double> distance_8 = {1};
    std::vector<double> label_8 = {0.0};

    // Setup parameters cryptocontext
    uint32_t multDepth = 39;
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


    vector<vector<Ciphertext<DCRTPoly>>> B = {c4,c3,c2,c1};
    vector<vector<Ciphertext<DCRTPoly>>> C = {c8,c7,c6,c5};

    for (int i=1; i<=6; i++){
        auto max = m_max(i,cc,B,C);
        Plaintext plaintext_dist;
        Plaintext plaintext_label;
        cc->Decrypt(keys.secretKey, max[0], &plaintext_dist);
        cc->Decrypt(keys.secretKey, max[1], &plaintext_label);
        cout << "Decrypted value (distance) of " << i << "-th max element: " << plaintext_dist << endl;
        cout << "Decrypted value (label) of " << i << "-th max element: " << plaintext_label << endl;
        }

    
    return 0;

}
    
    


    

    
   

    
