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

#include "ciphertext-ser.h"
#include "cryptocontext-ser.h"
#include "key/key-ser.h"
#include "scheme/ckksrns/ckksrns-ser.h"

using namespace lbcrypto;
using namespace std;


int main() {
    
    
    //STEP 0 - Read dataset
    //int num_train_samples = 16;

    ifstream file("C:/openfhe-development-main/src/pke/examples/heart_failure_clinical_records_dataset.csv");
        

        vector<vector<double>> dataset;

        string line;
        while (getline(file, line))
        {
            vector<double> row;
            stringstream ss(line);

            string cell;
            while (getline(ss, cell, ','))
            {
                    row.push_back(stod(cell));
                    
            }
            dataset.push_back(row);
        }

        file.close();
        dataset.assign(dataset.begin(), dataset.begin() + 16);

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


        //STEP 2 - Encriptar dades
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
        for (const auto& row : new_dataset) {
            vector<Ciphertext<DCRTPoly>> encrypted_row(row.size());
            for (size_t i = 0; i < row.size(); i++) {
                vector<double> values = { row[i][0] };
                Plaintext ptxt = cc->MakeCKKSPackedPlaintext(values);
                encrypted_row[i] = cc->Encrypt(keys.publicKey, ptxt);
            }
            encrypted_dataset.push_back(encrypted_row);
        }
        cout << "Dataset encrypted" << endl;

        cout << "Saving ciphertexts..." << endl;

        string DATAFOLDER = "C:/openfhe-development-main/src/pke/examples/encrypted_dataset/";
        for (size_t i = 0; i < encrypted_dataset.size(); i++) {
            for (size_t j = 0; j < encrypted_dataset[i].size(); j++) {
                string path = DATAFOLDER + "ciphertext" + to_string(i) + "_" + to_string(j) + ".txt";
                if (!Serial::SerializeToFile(path, encrypted_dataset[i][j], SerType::BINARY)) {
                    std::cerr << "Error writing serialization of ciphertext 3 to ciphertext3.txt" << std::endl;
                return 1;
                }
            }
        }
    
}

        
