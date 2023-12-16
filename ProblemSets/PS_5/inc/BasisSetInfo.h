//
// Created by Chongye Feng on 12/6/23.
//

#ifndef EHT_FINAL_BASISSETINFO_H
#define EHT_FINAL_BASISSETINFO_H

#include <vector>

struct BasisSetInfo {
    std::vector<double> exponents;
    std::vector<double> s_coefficients;
    std::vector<double> p_coefficients;
};

#endif //EHT_FINAL_BASISSETINFO_H
