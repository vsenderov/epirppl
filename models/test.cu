/* 
 *  models/CombineDS.cuh
 *
 *  Copyright (C) 2020-2021 Viktor Senderov, Joey Öhman, David Broman
 * 
 *
 *  CombineDS diversification model supports conditionally simulates
 *  several different types of evolution:
 *
 *    - cladogenetic (ClaDS-like) small changes in diversification
 *      rates, ClaDS versions 0-2.
 *	      
 *    - anagenetic small changes (happening on a single lineage)
 *
 *    - rare large shits (both anagenetic and cladogenetic)
 *
 *    - uncoupling of the turnover rate at the rare large shifts for
 *      ClaDS2
 *
 *
 *  This file needs to be included by a .cu file, containing the MAIN
 *  macro, needed global parameters, needed tuning parameters, and 
 *  the tree structure as a datatype.
 * 
 *  Needed global parameters:
 * 
 *    const floating_t k = 1;            // prior Γ-shape for λ
 *    const floating_t theta = 1;        // prior Γ-scale for λ
 *
 *    const floating_t kNu = 1;          // prior Γ-shape for ν
 *    const floating_t thetaNu = 0.5;    // prior Γ-shape for ν
 *
 *    const floating_t a_epsilon = 1;    // prior β-shape 1 for p_ε
 *    const floating_t b_epsilon = 100;  // prior β-shape 2 for p_ε
 *
 *    const floating_t m0 = 0;   // Hyper-param of prior for α and σ
 *    const floating_t v = 1;    // Hyper-param of prior for α and σ
 *    const floating_t a = 1.0;  // Hyper-param of prior for α and σ
 *    const floating_t b = 0.2;  // Hyper-param of prior for α and σ
 * 
 *  Needed tuning parameters:
 *
 *    #define M 20              // Number of subsamples to draw
 *    #define RARE_SHIFT false  // Activate rare shifts
 *    #define CLADS true        // Cladogenetic changes
 *    #define ANADS true        // Anagenetic changes
 *    #define UNCOUPLE true     // Uncouples turnover rate at rare shifts
 *    #define CLADS1 false      // ClaDS version: 0, 1, or 2, TODO 0
 *
 *  Tree selection, 3 steps:
 *
 *    #include "trees/cetaceans.cuh"       // (1)
 *    typedef cetaceans_87_tree_t tree_t;  // (2)
 *    const floating_t rhoConst = 1.00;    // (3) sampling rate
 *
 *  models/CombineDS.cuh defines the following BBLOCKS that can be included
 *  in the MAIN macro:
 *
 *    - simCombinedDS         (required)
 *
 *    - simTree               (required)
 *
 *    - conditionOnDetection  (optional, corrects for survivorship bias)
 *
 *    - sampleFinalLambda     (optional, samples the global parameters,
 *                             which have been delayed)
 *
 *    - saveResults           (optional callback, needs to be used in 
 *                             conjunction with sampleFinalLambda)
 */

/* Preamble */
#include <iostream>
#include <cstring>
#include <cassert>
#include <string>
#include <fstream>
#include <algorithm>
#include <random>

#include "inference/smc/smc.cuh"
#include "utils/math.cuh"
#include "utils/stack.cuh"
#include "dists/delayed.cuh"



/////////////////////////////////////////////////////////////////////////////



/* Program state */
struct progState_t {
  int x;
};


#define NUM_BBLOCKS 1

INIT_MODEL(progState_t)


/*
 * simCombineDS - required BBLOCK
 */
BBLOCK(test,
{
  //PSTATE.x = SAMPLE(binomial, 0.5, 100);
  for (int i = 0; i <= 6; i++)
    printf("%f\n", binomialScore(0.5, 5, i));
  NEXT = NULL;   
})
  

CALLBACK(COMPUTE, {
    floating_t sum = 0;
    for (int j = 0; j < N; j++)
      sum+=PSTATES[j].x;

    floating_t xbar = sum/N;

    floating_t sumsq;
    for (int j = 0; j < N; j++) {
      sumsq+=(PSTATES[j].x - xbar)*(PSTATES[j].x - xbar);
    }

    floating_t variance = sumsq/N;
    
    printf("mean: %f variance: %f\n", xbar, variance);
  })


MAIN({
    FIRST_BBLOCK(test);
    //SMC(COMPUTE)
    SMC(NULL);
})
