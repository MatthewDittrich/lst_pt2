#ifndef PT2OBJECT_H
#define PT2OBJECT_H

#include "rootReader.h"

class pt2Object {
public:
    int ls_idx  = -1;
    int pls_idx = -1;

    float ls_pt  = 0.0f;
    float ls_eta = 0.0f;
    float ls_phi = 0.0f;
    int ls_isFake = 1;
//    bool ls_isUsed = true;

    float pls_pt  = 0.0f;
    float pls_eta = 0.0f;
    float pls_phi = 0.0f;
    int pls_isFake = true;
//    bool pls_isUsed = true;

    float delta_pt = 0.0f;
    float delta_eta = 0.0f;
    float delta_phi = 0.0f;

//    float helix_dist1   = 999.99f;
//    float helix_dist2   = 999.99f;

    pt2Object() = default;

    pt2Object(const rootReader& reader, int lsIdx, int plsIdx);

private:
    void computeBasics();
//    void computeHelix();
};

#endif
