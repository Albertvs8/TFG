#define PROFILE

#include "openfhe.h"
#include <cmath>

using namespace lbcrypto;

//Plain inputs
double binomialCoefficient(int n, int k) {
    int result = 1;
    for (int i = 1; i <= k; i++) {
        result = result * (n - k + i) / i;
    }
    return result;
}

double f_n(int n, double x) {
    double sum = 0;
    for (int i = 0; i <= n; i++) {
        sum += (1.0 / pow(4, i)) * binomialCoefficient(2 * i, i) * x * pow(1 - pow(x, 2), i);
    }
    return sum;
}


double f_4(double x) {
    return (35.0 / 128.0) * pow(x, 9) - (180.0 / 128.0) * pow(x, 7) + (378.0 / 128.0) * pow(x, 5) - (420.0 / 128.0) * pow(x, 3) + (315.0 / 128.0) * x;
}

double g_4(double x) {
    return (46623.0 / 1024.0) * pow(x, 9) - (113492.0 / 1024.0) * pow(x, 7) + (97015.0 / 1024.0) * pow(x, 5) - (34974.0 / 1024.0) * pow(x, 3) + (5850.0 / 1024.0) * x;
}


double comparison_f(double a,double b,int n,int d) {
  double x = a - b;
  for (int i=1; i<=d; i++) {
    x = f_n(n, x);
  }
  return (1.0/2.0) * (x + 1.0);
}


double comparison_g(double a, double b, int n, int df, int dg) {
    double x = a - b;
    for (int i=1; i<=dg; i++) {
        x = g_4(x);
    }
    for (int i=1; i<=df; i++) {
        x = f_n(n, x);
    }
    return (1.0/2.0) * (x + 1.0);
}

//Encrypted inputs
std::vector<double> getInputVector(const std::string& name) {
    std::vector<double> vec;
    double input;

    std::cout << "Enter values for vector " << name << " (separated by commas): ";
    while (std::cin >> input) {
        vec.push_back(input);
        if (std::cin.peek() == ',') {
            std::cin.ignore();
        }
        else {
            break;
        }
    }
    return vec;
}

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




Ciphertext<DCRTPoly> homomorphic_comparison_f(const CryptoContext<DCRTPoly> &cc, const Ciphertext<DCRTPoly> &c1, const Ciphertext<DCRTPoly> &c2, int d, int vec_size) {
    auto x = cc->EvalSub(c1, c2);
    for (int i=1; i<=d; i++){
        x = f_4_homomorphic(cc,x,vec_size);
    }
    auto homomorphic_1 = cc->MakeCKKSPackedPlaintext(std::vector<double>(vec_size,{1.0}));
    auto homomorphic_05 = cc->MakeCKKSPackedPlaintext(std::vector<double>(vec_size,{0.5}));
    auto x_plus_1 = cc->EvalAdd(x, homomorphic_1);
    auto result = cc->EvalMult(x_plus_1, homomorphic_05);
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


int main() {
    //PLAIN INPUTS
    std::cout << "PLAIN INPUTS \n";

    double a;
    double b;
    int n=4;
    int d=9;
    int df=2;
    int dg=4;

    std::cout << "Enter a value for a: ";
    std::cin >> a;
    std::cout << "Enter a value for b: ";
    std::cin >> b;
    std::cout << ""<< std::endl;

    //Comparison using f
    auto start = std::chrono::high_resolution_clock::now();
    double result_f = comparison_f(a,b,n,d);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "Using only f: comp(" << a << ", " << b << ") = " << result_f << std::endl;
    std::cout << "Time taken: " << duration.count() << " microseconds" << std::endl;
    std::cout << ""<< std::endl;

    //Comparison using f and g
    start = std::chrono::high_resolution_clock::now();
    double result_g = comparison_g(a,b,n,df,dg);
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "Using g polynomials: comp(" << a << ", " << b << ") = " << result_g << std::endl;
    std::cout << "Time taken: " << duration.count() << " microseconds" << std::endl;
    std::cout << ""<< std::endl;



    //ENCRYPTED INPUTS
    std::cout << "ENCRYPTED INPUTS \n";

    // Read vectors
    std::vector<double> vector_a = getInputVector("a");
    std::vector<double> vector_b = getInputVector("b");
    std::cout << ""<< std::endl;

    // Error if vectors have different sizes
    if(vector_a.size() != vector_b.size()) {
        throw std::invalid_argument("Vectors are of different sizes");
    }
    
    int vec_size = vector_a.size();

    // Setup parameters cryptocontext
    uint32_t multDepth = 55;
    uint32_t scaleModSize = 50;
    uint32_t batchSize = vec_size;
    
    CCParams<CryptoContextCKKSRNS> parameters;
    
    parameters.SetMultiplicativeDepth(multDepth);
    parameters.SetScalingModSize(scaleModSize);
    parameters.SetScalingTechnique(FLEXIBLEAUTO);
    parameters.SetBatchSize(batchSize);

    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    // cc->Enable(ADVANCEDSHE); //for bootstrapping 
    // cc->Enable(FHE);   //for bootstrapping
    usint ringDim = cc->GetRingDimension();
    
    std::cout << "CKKS scheme is using ring dimension " << ringDim << std::endl;
    std::cout << ""<< std::endl;

    /* BOOTSTRAP SETUP
    std::vector<uint32_t> levelBudget = {4, 4};
    std::vector<uint32_t> bsgsDim = {0, 0};
    uint32_t numSlots = 8;
    cc->EvalBootstrapSetup(levelBudget,bsgsDim,numSlots);
    */

    auto keys = cc->KeyGen();
    cc->EvalMultKeyGen(keys.secretKey);

    // cc->EvalBootstrapKeyGen(keys.secretKey, 8);  // BOOTSTRAP keygen

    Plaintext ptxta = cc->MakeCKKSPackedPlaintext(vector_a);
    Plaintext ptxtb = cc->MakeCKKSPackedPlaintext(vector_b);

    std::cout << "Input a: " << ptxta << std::endl;
    std::cout << "Input b: " << ptxtb << std::endl;

    // Encrypt
    auto c1 = cc->Encrypt(keys.publicKey, ptxta);
    auto c2 = cc->Encrypt(keys.publicKey, ptxtb);

    Plaintext result_decrypted_f, result_decrypted_g;


    // Comparison using only f
    start = std::chrono::high_resolution_clock::now();
    auto result_f_enc = homomorphic_comparison_f(cc,c1,c2,d,vec_size);
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    cc->Decrypt(keys.secretKey, result_f_enc, &result_decrypted_f);
    result_decrypted_f->SetLength(batchSize);
    std::cout << "Using only f: comp( a ,  b ) = " << result_decrypted_f << std::endl;
    std::cout << "Time taken: " << std::chrono::duration<double>(duration).count() << " seconds" << std::endl;
    std::cout << ""<< std::endl;

    // Comparison using f and g
    start = std::chrono::high_resolution_clock::now();
    auto result_g_enc = homomorphic_comparison_g(cc,c1,c2,df,dg,vec_size);
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    cc->Decrypt(keys.secretKey, result_g_enc, &result_decrypted_g);
    result_decrypted_g->SetLength(batchSize);
    std::cout << "Using g polynomials: comp( a ,  b ) = " << result_decrypted_g << std::endl;
    std::cout << "Time taken: " << std::chrono::duration<double>(duration).count() << " seconds" << std::endl;

}
    
    


    

    
   

    
