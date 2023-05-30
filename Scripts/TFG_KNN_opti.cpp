#define PROFILE


#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>
#include <vector>
#include <string>
#include "openfhe.h"
#include <omp.h>
#include <cmath>
#include <stdexcept>

using namespace lbcrypto;
using namespace std;


//5.2.1
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
            double scaled_value;
            if (max_value-min_value==0){
                scaled_value = 0.0;
            } else{
                scaled_value = (0.288*(value - min_value)) / (max_value - min_value);
            }
            scaled_row[i][0] = scaled_value;
        }
        scaled_row.back()[0] = row.back()[0];
        scaled_dataset.push_back(scaled_row);
    }

    return scaled_dataset;
}



//5.2.3
// Function to compute the squared euclidean distance between two ciphertext arrays
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


//5.2.4
//f_4_homomorphic, f_3_homomorphic and f_2_homomorphic 
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
    auto result2 = cc->EvalAdd(term3, term4);
    result = cc->EvalAdd(result, result2);

    return result;
}


Ciphertext<DCRTPoly> f_2_homomorphic(const CryptoContext<DCRTPoly> &cc, const Ciphertext<DCRTPoly> &c, int vec_size) {
    auto c2 = cc->EvalMult(c, c); // c^2
    auto c3 = cc->EvalMult(c2, c); // c^3
    auto c5 = cc->EvalMult(c2, c3); // c^5

    auto term1 = cc->EvalMult(3.0/8.0, c5);
    auto term2 = cc->EvalMult(-10.0/8.0, c3);
    auto term3 = cc->EvalMult(15.0/8.0, c);

    auto result = cc->EvalAdd(term1, term2);
    result = cc->EvalAdd(result, term3);

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

    auto term1 = cc->EvalMult(-12860.0/1024.0, c7);
    auto term2 = cc->EvalMult(25614.0/1024.0, c5);
    auto term3 = cc->EvalMult(-16577.0/1024.0, c3);
    auto term4 = cc->EvalMult(4589.0/1024.0, c);

    auto result = cc->EvalAdd(term1, term2);
    auto result2 = cc->EvalAdd(term3, term4);
    result = cc->EvalAdd(result, result2);

    return result;
}

Ciphertext<DCRTPoly> g_2_homomorphic(const CryptoContext<DCRTPoly> &cc, const Ciphertext<DCRTPoly> &c, int vec_size) {
    auto c2 = cc->EvalMult(c, c); // c^2
    auto c3 = cc->EvalMult(c2, c); // c^3
    auto c5 = cc->EvalMult(c2, c3); // c^5

    auto term1 = cc->EvalMult(3796.0/1024.0, c5);
    auto term2 = cc->EvalMult(-6108.0/1024.0, c3);
    auto term3 = cc->EvalMult(3334.0/1024.0, c);

    auto result = cc->EvalAdd(term1, term2);
    result = cc->EvalAdd(result, term3);

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

    auto x_plus_1 = cc->EvalAdd(x, 1);
    auto result = cc->EvalMult(x_plus_1, 0.5);
    return result;
}


vector<vector<Ciphertext<DCRTPoly>>> compute_comparisons_parallel(const CryptoContext<DCRTPoly> &cc, const vector<vector<Ciphertext<DCRTPoly>>> &A, const Ciphertext<DCRTPoly> &c_enc_05) {
    //Number of compositions
    int df=2;
    int dg=3;
    int vec_size=1;

    vector<vector<Ciphertext<DCRTPoly>>> comparisonSetA(A.size(), vector<Ciphertext<DCRTPoly>>(A.size()));

    #pragma omp parallel for
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

    auto invComp = cc->EvalAdd(1,-comparison);

    auto first_term_dist = cc->EvalMult(comparison,F[0]);
    auto second_term_dist = cc->EvalMult(invComp,G[0]);
    auto result_dist = cc->EvalAdd(first_term_dist,second_term_dist);

    //L function calculation for label
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
    //int k = floor((s+t)/2);

    
    vector<vector<Ciphertext<DCRTPoly>>> Z(s+t);
    #pragma omp parallel for schedule(dynamic)
    for (int i=0; i<(s+t-2); i++){
        vector<Ciphertext<DCRTPoly>> zi = m_max(i+1,cc,B,C,comparisonSetBC);
        Z[i] = zi;
        
    }
    Z[s+t-1] = m_max(s+t,cc,B,C,comparisonSetBC);

    // sum all B
    Ciphertext<DCRTPoly> zi_dist=c0;
    Ciphertext<DCRTPoly> zi_label=c0;
    for (int i=0; i<s; i++){
        zi_dist = cc->EvalAdd(zi_dist,B[i][0]);
        zi_label = cc->EvalAdd(zi_label,B[i][1]);
    }
    // sum all C
    for (int i=0; i<t; i++){
        zi_dist = cc->EvalAdd(zi_dist,C[i][0]);
        zi_label = cc->EvalAdd(zi_label,C[i][1]);
    }

    //Subtract all Z (except last element)
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

    vector<vector<Ciphertext<DCRTPoly>>> B, C;

    #pragma omp task shared(B)
        {  
            B = merge_sort(cc,c0,vector<vector<Ciphertext<DCRTPoly>>>(A.begin(), A.begin()+s), comparisonSetB);
        }

    #pragma omp task shared(C)
        {
            C = merge_sort(cc,c0,vector<vector<Ciphertext<DCRTPoly>>>(A.begin()+s, A.end()), comparisonSetC);
        } 

    int k=A.size();

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

    auto return_value = merge(cc,c0,B,C,BCompC);
    return return_value;

}

int main() {
    
    
    //STEP 0 - Read dataset
    int num_train_samples = 15;
    int num_test_samples = 1;
    int max_0_test = floor(num_test_samples/2.0);
    int max_1_test = ceil(num_test_samples/2.0);

    for (int train_samples = num_train_samples; train_samples <= num_train_samples; train_samples++) {
        ifstream file("C:/openfhe-development-main/src/pke/examples/heart_failure_clinical_records_dataset.csv");
        
        int max_0 = floor(num_train_samples/2.0);
        int max_1 = ceil(num_train_samples/2.0);

        cout << "Train samples: " << num_train_samples << endl;
        cout << "Test samples: " << num_test_samples << endl;
        cout << endl;

        vector<vector<double>> dataset;
        vector<vector<double>> additional_records;
        int count_0 = 0, count_1 = 0; // counters for the number of records added
        int last_class = -1; // last class added to the dataset

        string line;
        while (getline(file, line))
        {
            vector<double> row;
            stringstream ss(line);

            string cell;
            int count = 0; // counter for the number of values pushed into the row vector
            while (getline(ss, cell, ','))
            {
                count++;
                if (count > -1) { // only push the last 4 values
                    row.push_back(stod(cell));
                }
            }

            // check if the last variable is 0 or 1 and add to the corresponding counter
            if (row.back() == 0 && count_0 < max_0) {
                if (last_class != 0) { // alternate classes
                    dataset.push_back(row);
                    count_0++;
                    last_class = 0;
                }
            }
            else if (row.back() == 1 && count_1 < max_1) {
                if (last_class != 1) { // alternate classes
                    dataset.push_back(row);
                    count_1++;
                    last_class = 1;
                }
            }
            else {
                additional_records.push_back(row);
            }
        }
        file.close();

        // add 2 additional records from one class and 2 from the other class
        int count_additional_0 = 0, count_additional_1 = 0;
        for (const auto& row : additional_records) {
            if (row.back() == 0 && count_additional_0 < max_0_test) {
                dataset.push_back(row);
                count_additional_0++;
            }
            else if (row.back() == 1 && count_additional_1 < max_1_test) {
                dataset.push_back(row);
                count_additional_1++;
            }
            if (count_additional_0 == max_0_test && count_additional_1 == max_1_test) {
                break; 
            }
        }

        // Transform each record from a vector to individual elements for encryption purposes 
        vector<vector<vector<double>>> new_dataset;
        transform(dataset.begin(), dataset.end(), back_inserter(new_dataset), 
                [](const vector<double>& row) {
                    vector<vector<double>> row_vectors;
                    transform(row.begin(), row.end(), back_inserter(row_vectors), 
                                [](double x) { return vector<double>{x}; });
                    return row_vectors;
                });
        
        // Print dataset
        cout << "Original dataset:" << endl;
        for (const auto& row : new_dataset) {
            for (const auto& element : row) {
                cout << "{" << element[0] << "}, ";
            }
            cout << endl;
        }
        cout << endl;

        
        //STEP 1 - Apply dataset scaling
        auto scaled_dataset = min_max_scaling(new_dataset);

        // Print scaled dataset
        cout << "Scaled dataset:" << endl;
        for (const auto& row : scaled_dataset) {
            for (const auto& element : row) {
                cout << "{" << element[0] << "}, ";
            }
            cout << endl;
        }
        cout << endl;


        //STEP 2 - Encriptar dades
        auto encrypt_start = std::chrono::high_resolution_clock::now();
        cout << "Encrypting dataset..." << endl;

        //2.1 Setup Parameters FHE
        uint32_t multDepth = 50;
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

        
        //Encryption of 0.5.  Variable auxiliar que es fa servir al pas 4.1
        std::vector<double> enc_05 = {0.5};
        Plaintext p_05 = cc->MakeCKKSPackedPlaintext(enc_05);
        auto c_enc_05 = cc->Encrypt(keys.publicKey, p_05);


        //Encryption of 0. Variable auxiliar que es fa servir al pas 4.2
        std::vector<double> enc_0 = {0};
        Plaintext p_0 = cc->MakeCKKSPackedPlaintext(enc_0);
        auto c_enc_0 = cc->Encrypt(keys.publicKey, p_0);

        auto encrypt_end = std::chrono::high_resolution_clock::now();
        auto duration_encrypt = std::chrono::duration_cast<std::chrono::microseconds>(encrypt_end - encrypt_start);


        std::cout << "Time taken in encryption process: " << std::chrono::duration<double>(duration_encrypt).count() << " seconds" << std::endl;
        cout << endl;

        // Separar train i test
        vector<vector<Ciphertext<DCRTPoly>>> train_dataset(encrypted_dataset.begin(), encrypted_dataset.begin() + num_train_samples);
        vector<vector<Ciphertext<DCRTPoly>>> test_dataset(encrypted_dataset.end() - num_test_samples, encrypted_dataset.end());

        for (size_t test_index=0; test_index<test_dataset.size(); test_index++){
            auto global_start = std::chrono::high_resolution_clock::now();
            //Agafar els features, no la label
            auto row_test = test_dataset[test_index];
            vector<Ciphertext<DCRTPoly>> features_test(row_test.begin(), row_test.end() - 1); 

            //Variable (array de distancies i labels) que s'ha d'ordenar
            vector<vector<Ciphertext<DCRTPoly>>> A;


            //STEP 3 - Calcular euclidean distance 
            cout << "Computing distances..." << endl;
            auto distance_start = std::chrono::high_resolution_clock::now();
            for (size_t train_index=0; train_index<train_dataset.size(); train_index++){
                auto row_train = train_dataset[train_index];
                vector<Ciphertext<DCRTPoly>> features_train(row_train.begin(), row_train.end() - 1 );
                Ciphertext<DCRTPoly> label_train = row_train[row_train.size()-1];


                Ciphertext<DCRTPoly> dist = euclidean_distance(cc,features_train,features_test);

                vector<Ciphertext<DCRTPoly>> dist_and_label = {dist, label_train};
                A.push_back(dist_and_label);
            }
            auto distance_stop = std::chrono::high_resolution_clock::now();
            auto distance_duration = std::chrono::duration_cast<std::chrono::microseconds>(distance_stop - distance_start);
            std::cout << "Time taken in computing distances: " << std::chrono::duration<double>(distance_duration).count() << " seconds" << std::endl;
            cout << endl;

            //STEP 4 - Ordenar euclidean distance (i respectiva label), fent servir merge_sort
            //4.1 Calcular pairwise comparisons
            cout << "Computing comparisons..." << endl;

            auto start = std::chrono::high_resolution_clock::now();
            vector<vector<Ciphertext<DCRTPoly>>> compSetA = compute_comparisons_parallel(cc,A,c_enc_05);
            auto stop_comp = std::chrono::high_resolution_clock::now();
            auto duration_comp = std::chrono::duration_cast<std::chrono::microseconds>(stop_comp - start);

            std::cout << "Time taken in computing comparisons: " << std::chrono::duration<double>(duration_comp).count() << " seconds" << std::endl;
            cout << endl;

            //4.2 Ordenar (Merge_sort)
            cout << "Ordering values..." << endl;

            vector<vector<Ciphertext<DCRTPoly>>> sorted_A = merge_sort(cc,c_enc_0,A,compSetA);
            auto stop_merge = std::chrono::high_resolution_clock::now();
            auto duration_merge = std::chrono::duration_cast<std::chrono::microseconds>(stop_merge - stop_comp);

            std::cout << "Time taken ordering: " << std::chrono::duration<double>(duration_merge).count() << " seconds" << std::endl;
            cout << endl;



            //STEP 5 - Agafar labels de les k distàncies més petites
            int k=3;
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
            double decrypt_sum_labels = result->GetRealPackedValue()[0];

            int predicted_label = -1;
            if (decrypt_sum_labels>k/2.0){
                predicted_label = 1;
            } else{
                predicted_label = 0;
            }


            cout << "Sum of labels: " << result << endl;
            cout << "Predicted label: " << predicted_label << endl;
             cout << endl;


            auto global_stop = std::chrono::high_resolution_clock::now();
            auto total_duration = std::chrono::duration_cast<std::chrono::microseconds>(global_stop - global_start);
            std::cout << "Total time taken in classifying: " << std::chrono::duration<double>(total_duration).count() << " seconds" << std::endl;

            
        }

                cout << "________________________________" << endl;
    }
    return 0;
}