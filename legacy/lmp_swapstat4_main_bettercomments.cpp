#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <list>
#include <string>
#include <fstream>
#include <algorithm>
#include <map>
#include <cmath>

#include "../../../general/stringext.h"
#include "../../../../aside/inih/cpp/INIReader.h"
#include "../../../math/basic_math/basic_math.h"
#include "../../../math/geometry/planar/vector.h"
#include "../../../math/geometry/planar/rect.h"
#include "../../../math/geometry/planar/circle.h"
#include "../../../math/stat/avg.h"
#include "../../../math/stat/avgerr.h"
#include "../../../math/stat/histogramm.h"
#include "../../../math/stat/histogramm2d.h"
#include "../../../math/stat/whist.h"
#include "../../../math/stat/whist2d.h"
#include "../../../math/stat/chist.h"
#include "../../../phys/atomic/pathanalysis.h"

using namespace std;

/**
* LAMMPS simulations parsed into the set of records with this structure
* Each swaprecord is an individual bistable bond snap event in simulations.
**/
struct swaprecord
{
    int nswap;
    int swapDirection;
    int atom1;
    int atom2;
    v3_math::spatial::vector<double> bondCenter;
    double timestep;
    double time;
    double xsize;
    double ysize;
    double surfaceArea;
    double _likelihood_ijdtheta;

    //int nsubavalanche;
    //int navalanche;
    //int nswapDirected;
    //double dtime;
    //double timeRatio;
    //double mds;
    //double mdsNormalized;
    //double distanceToPrevSwap;
    //double distanceToPrevAssocSwap;

    int set(const vector<string>& ss)
    {
        if (ss.size() != 24) return 1;
        nswap = v3_general::stoi(ss[0]);
        //nswapDirected = v3_general::stoi(ss[1]);
        timestep = v3_general::stod(ss[2]);
        time = v3_general::stod(ss[3]);
        //dtime = v3_general::stod(ss[4]);
        //timeRatio = v3_general::stod(ss[5]);
        swapDirection = v3_general::stoi(ss[6]);
        xsize = v3_general::stod(ss[7]);
        ysize = v3_general::stod(ss[8]);
        surfaceArea = v3_general::stod(ss[9]);
        //mds = v3_general::stod(ss[10]);
        //mdsNormalized = v3_general::stod(ss[11]);
        //nsubavalanche = v3_general::stoi(ss[14]);
        //navalanche = -1;
        bondCenter.x = v3_general::stod(ss[17]);
        bondCenter.y = v3_general::stod(ss[18]);
        bondCenter.z = v3_general::stod(ss[19]);
        //distanceToPrevSwap = v3_general::stoi(ss[20]);
        //distanceToPrevAssocSwap = v3_general::stoi(ss[21]);
        atom1 = v3_general::stoi(ss[22]);
        atom2 = v3_general::stoi(ss[23]);
        return 0;
    }
}; // struct swaprecord

/**
* A group of individual bistable bond snaps happening close in space and time
* is called subavalanche or cascade.
**/
struct subavalanche
{
    int subavalancheID;
    vector<pair<double, int>> swaprecords_time_ID;

    int nswaps;
    double subavTimeP50;
    double subavTimeStart;
    double subavTimeEnd;
    double subavDuration;
    double subavCenterX;
    double subavCenterY;
    double subavRg;
    double subavI1, subavI2;
    double subavDiameter;
};

/** 
* A group of subsequent cascades represents thermal avalanche, called here just
* avalanche.
**/
struct avalanche
{
    int avalancheID;
    vector<pair<double, int>> subavalanches_timeP50_ID;
    vector<pair<double, int>> subavalanches_timeStart_ID;
    vector<pair<double, int>> swaprecords_time_ID;

    int avNsubavs;
    int avNswaps;
    int avNswapsDir;
    double avTimeStart;
    double avTimeEnd;
    double avTimeAvg;
    double avThetaAvg;
    double avTimeP50;
    double avDuration;
    double avThetaDuration;
    double avTg;
    double avThetag;
    double avRg;
    double avDiameter;
};


int main(int argc, char* argv[])
{

/**
* Technical block for initialization and data parsing
**/

    /** $ Read command line arguments **/
    if (argc < 3)
    {
        cout << "\nNot enough arguments!";
        return 1;
    }

    string swapfnfile = argv[1];
    string paramfname = argv[2];
    string outfname = argv[3];
    /** / Read command line arguments **/


    /** $ Read ini file **/
    cout << "\nRead ini file...";
    // Open file
    INIReader ini(paramfname);

    // Check parsing results
    switch (ini.ParseError())
    {
    case 0: // parsed successfully
        break;
    default:
        cout << "Error parsing parameters file \"" << paramfname << "\"\n";
        return 1;
    }

    // System
    double system_posLeft_x = ini.GetReal("System", "system_posLeft_x", 0.00);
    double system_posLeft_y = ini.GetReal("System", "system_posLeft_y", 0.00);

    // SubavalancheDetection
    double avgBondLength = ini.GetReal("SubavalancheDetection", "avgBondLength", 0.00);
    double subavalanche_dt_threshold = ini.GetReal("SubavalancheDetection", "subavalanche_dt_threshold", 10.00);
    double subavalanche_dr_threshold = ini.GetReal("SubavalancheDetection", "subavalanche_dr_threshold", 100.00);

    // Avalanche detection
    // Avalanche detection mechanism:
    // 0 - by dt/t cutoff
    // 1 - by dt/t and r cutoffs, closest in space wins
    // 2 - by dt/t and and r cutoffs, bigger wins
    // 3 - by dt/t and and r cutoffs plus minimum dt of the avalanche (to prevent small dt/t cutoffs in the beginning)
    // 4 - by dt/t and r/dt > v cutoffs, bigger wins
    // 5 - by dt/t and r/dt > v cutoffs, biggest choice
    int avalancheDetectionMechanism = ini.GetInteger("AvalancheDetection", "avalancheDetectionMechanism", 0); 
    int avDetection_subavNswaps_threshold = ini.GetInteger("AvalancheDetection", "avDetection_subavNswaps_threshold", 1); // subavalanche nswaps threshold for mechanism == 2
    double avDetection_dtt_threshold = ini.GetReal("AvalancheDetection", "avDetection_dtt_threshold", 0.01); // dt/t threshold for mechanism == 0 || 1 || 2 || 3 || 4 || 5
    double avDetection_dt_minimum = ini.GetReal("AvalancheDetection", "avDetection_dt_minimum", 0.01); // dt minimum for mechanism == 3
    double avDetection_dr_threshold = ini.GetReal("AvalancheDetection", "avDetection_dr_threshold", 100.00); // r threshold for mechanism == 1 || 2 || 3
    double avDetection_velocity_threshold = ini.GetReal("AvalancheDetection", "avDetection_velocity_threshold", 20.00); // velocity threshold for mechanism ==  4 || 5

    // Statistics
    double analysis_time_min = ini.GetReal("Statistics", "analysis_time_min", 0.00);
    double analysis_time_max = ini.GetReal("Statistics", "analysis_time_max", -1.00);

    // Swaps statistics
    double swapStat_epochStart = ini.GetReal("SwapStat", "epochStart", 10.00);
    double swapStat_epochStart_log10 = log10(swapStat_epochStart);
    int swapStat_nepoch = ini.GetInteger("SwapStat", "nepoch", 6);
    int swapStat_nepochfrac = ini.GetInteger("SwapStat", "nepochfrac", 2);
    double swapStat_dtdr_logdt_min = ini.GetReal("SwapStat", "dtdr_logdt_min", -4.0);
    double swapStat_dtdr_logdt_max = ini.GetReal("SwapStat", "dtdr_logdt_max", 6.0);
    int swapStat_dtdr_logdt_nbins = ini.GetInteger("SwapStat", "dtdr_logdt_nbins", 20);
    double swapStat_dtdr_logdr_min = ini.GetReal("SwapStat", "dtdr_logdr_min", 0.0);
    double swapStat_dtdr_logdr_max = ini.GetReal("SwapStat", "dtdr_logdr_max", 5.0);
    int swapStat_dtdr_logdr_nbins = ini.GetInteger("SwapStat", "dtdr_logdr_nbins", 20);
    double swapStat_dr_min = ini.GetReal("SwapStat", "dr_min", 0.0);
    double swapStat_dr_max = ini.GetReal("SwapStat", "dr_max", 500.0);
    int swapStat_dr_nbins = ini.GetInteger("SwapStat", "dr_nbins", 20);
    double swapStat_logdr_min = ini.GetReal("SwapStat", "logdr_min", 1.0);
    double swapStat_logdr_max = ini.GetReal("SwapStat", "logdr_max", 500.0);
    int swapStat_logdr_nbins = ini.GetInteger("SwapStat", "logdr_nbins", 20);

    double swapStat_ij_dtdr_logdt_min = ini.GetReal("SwapStat", "ij_dtdr_logdt_min", -4.0);
    double swapStat_ij_dtdr_logdt_max = ini.GetReal("SwapStat", "ij_dtdr_logdt_max", 7.0);
    int swapStat_ij_dtdr_logdt_nbins = ini.GetInteger("SwapStat", "ij_dtdr_logdt_nbins", 22);
    double swapStat_ij_dtdr_logdr_min = ini.GetReal("SwapStat", "ij_dtdr_logdr_min", 0.0);
    double swapStat_ij_dtdr_logdr_max = ini.GetReal("SwapStat", "ij_dtdr_logdr_max", 5.0);
    int swapStat_ij_dtdr_logdr_nbins = ini.GetInteger("SwapStat", "ij_dtdr_logdr_nbins", 20);
    double swapStat_ij_dt_min = ini.GetReal("SwapStat", "ij_dt_min", 0.0);
    double swapStat_ij_dt_max = ini.GetReal("SwapStat", "ij_dt_max", 5.0);
    int swapStat_ij_dt_nbins = ini.GetInteger("SwapStat", "ij_dt_nbins", 20);
    double swapStat_ij_dtheta_min = ini.GetReal("SwapStat", "ij_dtheta_min", 0.0);
    double swapStat_ij_dtheta_max = ini.GetReal("SwapStat", "ij_dtheta_max", 5.0);
    int swapStat_ij_dtheta_nbins = ini.GetInteger("SwapStat", "ij_dtheta_nbins", 20);
    double swapStat_ij_dr_min = ini.GetReal("SwapStat", "ij_dr_min", 0.0);
    double swapStat_ij_dr_max = ini.GetReal("SwapStat", "ij_dr_max", 5.0);
    int swapStat_ij_dr_nbins = ini.GetInteger("SwapStat", "ij_dr_nbins", 20);
    
    // Subavalanche statistics
    int simultaneousSubavalanches_mode = ini.GetInteger("SubavalancheStat", "simultaneousSubavalanches_mode", 0);
    double subavStat_velocity = ini.GetReal("SubavalancheStat", "subavStat_velocity", 10.0);
    double subavStat_nswaps_log10_min = ini.GetReal("SubavalancheStat", "nswaps_log10_min", 0.00);
    double subavStat_nswaps_log10_max = ini.GetReal("SubavalancheStat", "nswaps_log10_max", 0.00);
    int subavStat_nswaps_log10_nbins = ini.GetInteger("SubavalancheStat", "nswaps_log10_nbins", 0.00);
    double subavStat_duration_log10_min = ini.GetReal("SubavalancheStat", "duration_log10_min", 0.00);
    double subavStat_duration_log10_max = ini.GetReal("SubavalancheStat", "duration_log10_max", 0.00);
    int subavStat_duration_log10_nbins = ini.GetInteger("SubavalancheStat", "duration_log10_nbins", 0.00);
    double subavStat_swaprate_log10_min = ini.GetReal("SubavalancheStat", "swaprate_log10_min", 0.00);
    double subavStat_swaprate_log10_max = ini.GetReal("SubavalancheStat", "swaprate_log10_max", 0.00);
    int subavStat_swaprate_log10_nbins = ini.GetInteger("SubavalancheStat", "swaprate_log10_nbins", 0.00);
    double subavStat_diameter_log10_min = ini.GetReal("SubavalancheStat", "diameter_log10_min", 0.00);
    double subavStat_diameter_log10_max = ini.GetReal("SubavalancheStat", "diameter_log10_max", 0.00);
    int subavStat_diameter_log10_nbins = ini.GetInteger("SubavalancheStat", "diameter_log10_nbins", 0.00);
    double subavStat_rg_log10_min = ini.GetReal("SubavalancheStat", "rg_log10_min", 0.00);
    double subavStat_rg_log10_max = ini.GetReal("SubavalancheStat", "rg_log10_max", 0.00);
    int subavStat_rg_log10_nbins = ini.GetInteger("SubavalancheStat", "rg_log10_nbins", 0.00);
    double subavStat_I12_log10_min = ini.GetReal("SubavalancheStat", "I12_log10_min", 0.00);
    double subavStat_I12_log10_max = ini.GetReal("SubavalancheStat", "I12_log10_max", 0.00);
    int subavStat_I12_log10_nbins = ini.GetInteger("SubavalancheStat", "I12_log10_nbins", 0.00);
    /*double density_log10_min = ini.GetReal("SubavalancheStat", "density_log10_min", 0.00);
    double density_log10_max = ini.GetReal("SubavalancheStat", "density_log10_max", 0.00);
    int density_log10_nbins = ini.GetInteger("SubavalancheStat", "density_log10_nbins", 0.00);*/
    double subavStat_dtP50_log10_min = ini.GetReal("SubavalancheStat", "dtP50_log10_min", 0.00);
    double subavStat_dtP50_log10_max = ini.GetReal("SubavalancheStat", "dtP50_log10_max", 0.00);
    int subavStat_dtP50_log10_nbins = ini.GetInteger("SubavalancheStat", "dtP50_log10_nbins", 0.00);
    double subavStat_dttP50_log10_min = ini.GetReal("SubavalancheStat", "dttP50_log10_min", 0.00);
    double subavStat_dttP50_log10_max = ini.GetReal("SubavalancheStat", "dttP50_log10_max", 0.00);
    int subavStat_dttP50_log10_nbins = ini.GetInteger("SubavalancheStat", "dttP50_log10_nbins", 0.00);
    double subavStat_dr_min = ini.GetReal("SubavalancheStat", "dr_min", 0.00);
    double subavStat_dr_max = ini.GetReal("SubavalancheStat", "dr_max", 0.00);
    int subavStat_dr_nbins = ini.GetInteger("SubavalancheStat", "dr_nbins", 0.00);
    double subavStat_dr_log10_min = ini.GetReal("SubavalancheStat", "dr_log10_min", 0.00);
    double subavStat_dr_log10_max = ini.GetReal("SubavalancheStat", "dr_log10_max", 0.00);
    int subavStat_dr_log10_nbins = ini.GetInteger("SubavalancheStat", "dr_log10_nbins", 0.00);
    double subavStat_ij_dtP50_log10_min = ini.GetReal("SubavalancheStat", "ij_dtP50_log10_min", 0.00);
    double subavStat_ij_dtP50_log10_max = ini.GetReal("SubavalancheStat", "ij_dtP50_log10_max", 0.00);
    int subavStat_ij_dtP50_log10_nbins = ini.GetInteger("SubavalancheStat", "ij_dtP50_log10_nbins", 0.00);
    double subavStat_ij_dttP50_log10_min = ini.GetReal("SubavalancheStat", "ij_dttP50_log10_min", 0.00);
    double subavStat_ij_dttP50_log10_max = ini.GetReal("SubavalancheStat", "ij_dttP50_log10_max", 0.00);
    int subavStat_ij_dttP50_log10_nbins = ini.GetInteger("SubavalancheStat", "ij_dttP50_log10_nbins", 0.00);
    double subavStat_ij_dr_min = ini.GetReal("SubavalancheStat", "ij_dr_min", 0.00);
    double subavStat_ij_dr_max = ini.GetReal("SubavalancheStat", "ij_dr_max", 0.00);
    int subavStat_ij_dr_nbins = ini.GetInteger("SubavalancheStat", "ij_dr_nbins", 0.00);
    v3_phys::path_analyser_setup subav_path_setup;
    // dr distribution
    subav_path_setup.dr_log = ini.GetBoolean("SubavalancheStat", "subavPath_dr_log", 0);
    subav_path_setup.dr_min = ini.GetReal("SubavalancheStat", "subavPath_dr_min", 1.0);
    subav_path_setup.dr_max = ini.GetReal("SubavalancheStat", "subavPath_dr_max", 2.0);
    subav_path_setup.dr_nbins = ini.GetInteger("SubavalancheStat", "subavPath_dr_nbins", 1);
    // dt distribution
    subav_path_setup.dt_log = ini.GetBoolean("SubavalancheStat", "subavPath_dt_log", 0);
    subav_path_setup.dt_min = ini.GetReal("SubavalancheStat", "subavPath_dt_min", 1.0);
    subav_path_setup.dt_max = ini.GetReal("SubavalancheStat", "subavPath_dt_max", 2.0);
    subav_path_setup.dt_nbins = ini.GetInteger("SubavalancheStat", "subavPath_dt_nbins", 1);
    // mean square dpath-tau
    subav_path_setup.dpath_tau_taulog = ini.GetBoolean("SubavalancheStat", "subavPath_dpath_tau_taulog", 0);
    subav_path_setup.dpath_tau_taumin = ini.GetReal("SubavalancheStat", "subavPath_dpath_tau_taumin", 1.0);
    subav_path_setup.dpath_tau_taumax = ini.GetReal("SubavalancheStat", "subavPath_dpath_tau_taumax", 2.0);
    subav_path_setup.dpath_tau_taunbins = ini.GetInteger("SubavalancheStat", "subavPath_dpath_tau_taunbins", 1);
    // mean square dr-tau
    subav_path_setup.dr_tau_taulog = ini.GetBoolean("SubavalancheStat", "subavPath_dr_tau_taulog", 0);
    subav_path_setup.dr_tau_taumin = ini.GetReal("SubavalancheStat", "subavPath_dr_tau_taumin", 1.0);
    subav_path_setup.dr_tau_taumax = ini.GetReal("SubavalancheStat", "subavPath_dr_tau_taumax", 2.0);
    subav_path_setup.dr_tau_taunbins = ini.GetInteger("SubavalancheStat", "subavPath_dr_tau_taunbins", 1);
    // mean square dr-dpath
    subav_path_setup.dr_dpath_dpathlog = ini.GetBoolean("SubavalancheStat", "subavPath_dr_dpath_dpathlog", 0);
    subav_path_setup.dr_dpath_dpathmin = ini.GetReal("SubavalancheStat", "subavPath_dr_dpath_dpathmin", 1.0);
    subav_path_setup.dr_dpath_dpathmax = ini.GetReal("SubavalancheStat", "subavPath_dr_dpath_dpathmax", 2.0);
    subav_path_setup.dr_dpath_dpathnbins = ini.GetInteger("SubavalancheStat", "subavPath_dr_dpath_dpathnbins", 1);
    // Fractal dimentions
    int subav_fractDim_nLengths = ini.GetInteger("SubavalancheStat", "subav_fractDim_nLengths", 0);
    int subav_fractDim_nswapsmin = ini.GetInteger("SubavalancheStat", "subav_fractDim_nswapsmin", 0);
    int subav_fractDim_nswapsmax = ini.GetInteger("SubavalancheStat", "subav_fractDim_nswapsmax", 10000);
    double subav_fractDim_Rgmin = ini.GetReal("SubavalancheStat", "subav_fractDim_Rgmin", 0.00);
    double subav_fractDim_Rgmax = ini.GetReal("SubavalancheStat", "subav_fractDim_Rgmax", 1000.00);
    double subav_fractDim_logS_min = ini.GetReal("SubavalancheStat", "subav_fractDim_logS_min", 0.00);
    double subav_fractDim_logS_max = ini.GetReal("SubavalancheStat", "subav_fractDim_logS_max", 8.00);
    int subav_fractDim_logS_nbins = ini.GetInteger("SubavalancheStat", "subav_fractDim_logS_nbins", 16);
    double subav_fractDim_logL_min = ini.GetReal("SubavalancheStat", "subav_fractDim_logL_min", 0.00);
    double subav_fractDim_logL_max = ini.GetReal("SubavalancheStat", "subav_fractDim_logL_max", 8.00);
    int subav_fractDim_logL_nbins = ini.GetInteger("SubavalancheStat", "subav_fractDim_logL_nbins", 16);
    double subav_fractDim_logStoL_min = ini.GetReal("SubavalancheStat", "subav_fractDim_logStoL_min", 0.00);
    double subav_fractDim_logStoL_max = ini.GetReal("SubavalancheStat", "subav_fractDim_logStoL_max", 8.00);
    int subav_fractDim_logStoL_nbins = ini.GetInteger("SubavalancheStat", "subav_fractDim_logStoL_nbins", 16);

    // Avalanche statistics
    // Single avalanche
    int avStat_nswaps_threshold = ini.GetInteger("AvalancheStat", "nswaps_threshold", 0);
    double avStat_nswaps_log10_min = ini.GetReal("AvalancheStat", "nswaps_log10_min", 0.00);
    double avStat_nswaps_log10_max = ini.GetReal("AvalancheStat", "nswaps_log10_max", 0.00);
    int avStat_nswaps_log10_nbins = ini.GetInteger("AvalancheStat", "nswaps_log10_nbins", 0.00);
    double avStat_nsubavs_min = ini.GetReal("AvalancheStat", "nsubavs_min", 1.0);
    double avStat_nsubavs_max = ini.GetReal("AvalancheStat", "nsubavs_max", 2.0);
    int avStat_nsubavs_nbins = ini.GetInteger("AvalancheStat", "nsubavs_nbins", 1);
    double avStat_duration_log10_min = ini.GetReal("AvalancheStat", "duration_log10_min", 1.0);
    double avStat_duration_log10_max = ini.GetReal("AvalancheStat", "duration_log10_max", 2.0);
    int avStat_duration_log10_nbins = ini.GetInteger("AvalancheStat", "duration_log10_nbins", 1);
    double avStat_tg_log10_min = ini.GetReal("AvalancheStat", "tg_log10_min", 1.0);
    double avStat_tg_log10_max = ini.GetReal("AvalancheStat", "tg_log10_max", 2.0);
    int avStat_tg_log10_nbins = ini.GetInteger("AvalancheStat", "tg_log10_nbins", 1);
    double avStat_dttg_log10_min = ini.GetReal("AvalancheStat", "dttg_log10_min", 1.0);
    double avStat_dttg_log10_max = ini.GetReal("AvalancheStat", "dttg_log10_max", 2.0);
    int avStat_dttg_log10_nbins = ini.GetInteger("AvalancheStat", "dttg_log10_nbins", 1);
    double avStat_Rg_log10_min = ini.GetReal("AvalancheStat", "Rg_log10_min", 1.0);
    double avStat_Rg_log10_max = ini.GetReal("AvalancheStat", "Rg_log10_max", 2.0);
    int avStat_Rg_log10_nbins = ini.GetInteger("AvalancheStat", "Rg_log10_nbins", 1);
    double avStat_dtheta_log10_min = ini.GetReal("AvalancheStat", "dtheta_log10_min", 1.0);
    double avStat_dtheta_log10_max = ini.GetReal("AvalancheStat", "dtheta_log10_max", 2.0);
    int avStat_dtheta_log10_nbins = ini.GetInteger("AvalancheStat", "dtheta_log10_nbins", 1);
    double avStat_ij_dtheta_log10_min = ini.GetReal("AvalancheStat", "ij_dtheta_log10_min", 0.00);
    double avStat_ij_dtheta_log10_max = ini.GetReal("AvalancheStat", "ij_dtheta_log10_max", 0.00);
    int avStat_ij_dtheta_log10_nbins = ini.GetInteger("AvalancheStat", "ij_dtheta_log10_nbins", 0.00);
    int avStructure_nswaps_threshold = ini.GetInteger("AvalancheStat", "avStructure_nswaps_threshold", 0);
    double avStructure_dt = ini.GetReal("AvalancheStat", "avStructure_dt", 1.0);
    double avStructure_dlogt = ini.GetReal("AvalancheStat", "avStructure_dlogt", 1.0);
    double avStructure_dlogtnorm = ini.GetReal("AvalancheStat", "avStructure_dlogtnorm", 0.1);
    int av_fractDim_nswapsmin = ini.GetInteger("AvalancheStat", "av_fractDim_nswapsmin", 0);
    int av_fractDim_nswapsmax = ini.GetInteger("AvalancheStat", "av_fractDim_nswapsmax", 10000);
    
    // Path
    v3_phys::path_analyser_setup av_path_setup;
    // dr distribution
    av_path_setup.dr_log = ini.GetBoolean("AvalancheStat", "avPath_dr_log", 0);
    av_path_setup.dr_min = ini.GetReal("AvalancheStat", "avPath_dr_min", 1.0);
    av_path_setup.dr_max = ini.GetReal("AvalancheStat", "avPath_dr_max", 2.0);
    av_path_setup.dr_nbins = ini.GetInteger("AvalancheStat", "avPath_dr_nbins", 1);
    // dt distribution
    av_path_setup.dt_log = ini.GetBoolean("AvalancheStat", "avPath_dt_log", 0);
    av_path_setup.dt_min = ini.GetReal("AvalancheStat", "avPath_dt_min", 1.0);
    av_path_setup.dt_max = ini.GetReal("AvalancheStat", "avPath_dt_max", 2.0);
    av_path_setup.dt_nbins = ini.GetInteger("AvalancheStat", "avPath_dt_nbins", 1);
    // mean square dpath-tau
    av_path_setup.dpath_tau_taulog = ini.GetBoolean("AvalancheStat", "avPath_dpath_tau_taulog", 0);
    av_path_setup.dpath_tau_taumin = ini.GetReal("AvalancheStat", "avPath_dpath_tau_taumin", 1.0);
    av_path_setup.dpath_tau_taumax = ini.GetReal("AvalancheStat", "avPath_dpath_tau_taumax", 2.0);
    av_path_setup.dpath_tau_taunbins = ini.GetInteger("AvalancheStat", "avPath_dpath_tau_taunbins", 1);
    // mean square dr-tau
    av_path_setup.dr_tau_taulog = ini.GetBoolean("AvalancheStat", "avPath_dr_tau_taulog", 0);
    av_path_setup.dr_tau_taumin = ini.GetReal("AvalancheStat", "avPath_dr_tau_taumin", 1.0);
    av_path_setup.dr_tau_taumax = ini.GetReal("AvalancheStat", "avPath_dr_tau_taumax", 2.0);
    av_path_setup.dr_tau_taunbins = ini.GetInteger("AvalancheStat", "avPath_dr_tau_taunbins", 1);
    // mean square dr-dpath
    av_path_setup.dr_dpath_dpathlog = ini.GetBoolean("AvalancheStat", "avPath_dr_dpath_dpathlog", 0);
    av_path_setup.dr_dpath_dpathmin = ini.GetReal("AvalancheStat", "avPath_dr_dpath_dpathmin", 1.0);
    av_path_setup.dr_dpath_dpathmax = ini.GetReal("AvalancheStat", "avPath_dr_dpath_dpathmax", 2.0);
    av_path_setup.dr_dpath_dpathnbins = ini.GetInteger("AvalancheStat", "avPath_dr_dpath_dpathnbins", 1);
    // Between avalanches
    double avStat_dtP50_log10_min = ini.GetReal("AvalancheStat", "dtP50_log10_min", 0.00);
    double avStat_dtP50_log10_max = ini.GetReal("AvalancheStat", "dtP50_log10_max", 0.00);
    int avStat_dtP50_log10_nbins = ini.GetInteger("AvalancheStat", "dtP50_log10_nbins", 0.00);
    double avStat_dttP50_log10_min = ini.GetReal("AvalancheStat", "dttP50_log10_min", 0.00);
    double avStat_dttP50_log10_max = ini.GetReal("AvalancheStat", "dttP50_log10_max", 0.00);
    int avStat_dttP50_log10_nbins = ini.GetInteger("AvalancheStat", "dttP50_log10_nbins", 0.00);
    double avStat_dthetaP50_log10_min = ini.GetReal("AvalancheStat", "dthetaP50_log10_min", 0.00);
    double avStat_dthetaP50_log10_max = ini.GetReal("AvalancheStat", "dthetaP50_log10_max", 0.00);
    int avStat_dthetaP50_log10_nbins = ini.GetInteger("AvalancheStat", "dthetaP50_log10_nbins", 0.00);
    
    cout << "done.\n";
    /** / Read ini file **/

    /** $ Read swap file names **/
    int nfiles;
    vector<string> swapfnames;
    vector<double> maxsimtime;
    ifstream swapfnf(swapfnfile.c_str());
    if (!swapfnf)
    {
        cout << "\nCan not open file \"" << swapfnfile << "\"";
        return 1;
    }
    string s;
    double mst;
    swapfnf >> s >> nfiles;
    for (int q = 0; q < nfiles; q++)
    {
        swapfnf >> s >> mst;
        swapfnames.push_back(s);
        maxsimtime.push_back(mst);
    }
    swapfnf.close();
    cout << "done.\n";
    /** / Read swap file names **/


    /** $ Read swap files **/
    cout << "Read swap files...\n";
    vector<vector<swaprecord>> swaprecords(nfiles);
    for (int nf = 0; nf < nfiles; nf++)
    {
        cout << "File " << swapfnames[nf] << "\n";
        ifstream inpf(swapfnames[nf].c_str());
        if (!inpf)
        {
            cout << "\nCan not open file \"" << swapfnames[nf] << "\"";
            continue;
        }
        int nlines = 0;
        int nerrors = 0;
        string line;
        vector<string> tokens;
        swaprecord r;
        while (getline(inpf, line))
        { // '\n' is the default delimiter
            nlines++;
            if (nlines == 1) continue;
            istringstream iss(line);
            tokens.clear();
            while (getline(iss, s, '\t'))   // but we can specify a different one
                tokens.push_back(s);
            nerrors += r.set(tokens);
            swaprecords[nf].push_back(r);
        }
        cout << "Records read: " << nlines - 1 << ", errors " << nerrors << "\n";
        inpf.close();
    }
    cout << "done.\n";
    /** / Read swap files **/


/**
* This block combines bond snaps into subavalanches (cascades)
**/
    /** $ Detect subavalanches **/
    cout << "Detect subavalanches...\n";
    vector<vector<int>> swap_subavalancheID(nfiles);
    vector<vector<subavalanche>> subavalanches(nfiles);
    vector<vector<int>> subavalanches_nsimultaneous(nfiles);

    if (true /*subavalancheDetectionMechanism == 0*/)
    {
        for (int nf = 0; nf < nfiles; nf++)
        {
            int nsubavalanches = 0;
            int nswaps = swaprecords[nf].size();
            swap_subavalancheID[nf].resize(nswaps);
            for (int q = 0; q < nswaps; q++) swap_subavalancheID[nf][q] = -1;

            vector<vector<int>> current_subavalanches; // swap IDs in currently open subavalanches
            for (int ns = 0; ns < nswaps; ns++)
            {
                //cout << "\nswap" << ns;
                //cout.flush();
                // Calculate dt and dr for each of the current subavalanches
                int nsubav = current_subavalanches.size();
                vector<double> dt_last(nsubav), dr_threshold_min(nsubav), dt_threshold_min(nsubav);
                for (int q = 0; q < nsubav; q++)
                {
                    int saSize = current_subavalanches[q].size();
                    dt_last[q] = swaprecords[nf][ns].time - swaprecords[nf][current_subavalanches[q][saSize - 1]].time;

                    dr_threshold_min[q] = -1.00;
                    for (int w = 0; w < saSize; w++)
                    {
                        double dt = swaprecords[nf][ns].time - swaprecords[nf][current_subavalanches[q][saSize - 1 - w]].time;
                        if (dt > subavalanche_dt_threshold) break;
                        double dr = (swaprecords[nf][ns].bondCenter - swaprecords[nf][current_subavalanches[q][saSize - 1 - w]].bondCenter).length();
                        if (dr_threshold_min[q] >= 0.00)
                        {
                            if (dr_threshold_min[q] > dr)
                            {
                                dr_threshold_min[q] = dr;
                                dt_threshold_min[q] = dt;
                            }
                        }
                        else
                        {
                            dr_threshold_min[q] = dr;
                            dt_threshold_min[q] = dt;
                        }
                    }
                }
                //cout << " dtdr" << nsubav;
                //cout.flush();

                // Determine the closest subavalanche
                int closestSubav = -1;
                double dr_min = 0.00;
                for (int q = 0; q < nsubav; q++)
                {
                    if (dr_threshold_min[q] >= 0.00)
                    {
                        if (dr_threshold_min[q] < subavalanche_dr_threshold)
                        {
                            if (closestSubav < 0)
                            {
                                dr_min = dr_threshold_min[q];
                                closestSubav = q;
                            }
                            else if (dr_min > dr_threshold_min[q])
                            {
                                dr_min = dr_threshold_min[q];
                                closestSubav = q;
                            }
                        }
                    }
                }
                //cout << " cls" << closestSubav;
                //cout.flush();

                // Check if the new event breaks any of the subavalanche sequences
                vector<bool> isSubavStopped(nsubav);
                int nstop = 0;
                for (int q = 0; q < nsubav; q++)
                {
                    isSubavStopped[q] = false;
                    if (dt_last[q] > subavalanche_dt_threshold)
                    {
                        isSubavStopped[q] = true;
                        nstop++;
                        continue;
                    }
                }
                //cout << " stop" << nstop;
                //cout.flush();

                // Add the event to the closest subavalanche (if exists)
                // or to the new avalanche
                if (closestSubav >= 0)
                {
                    current_subavalanches[closestSubav].push_back(ns);
                    //current_subavalanche_dtime.add(swap_dtime_sorted[ns]);
                    //swap_subavalanche_avgdt[ns] = current_subavalanche_dtime.getavg();
                }
                else
                {
                    vector<int> newsa(1);
                    newsa[0] = ns;
                    current_subavalanches.push_back(newsa);
                    isSubavStopped.push_back(false);
                }
                //cout << " add" << current_subavalanches.size() - nsubav;
                //cout.flush();

                // Store & forget stopped current subavalanches
                bool tryErasing = true;
                while (tryErasing)
                {
                    tryErasing = false;
                    for (int q = 0; q < current_subavalanches.size(); q++)
                    {
                        if (isSubavStopped[q])
                        {
                            //cout << "\nStore subav " << q << " of " << current_subavalanches.size()
                            //    << " with " << current_subavalanches[q].size() << " events";
                            //cout.flush();
                            subavalanche newsa;
                            newsa.subavalancheID = nsubavalanches;
                            for (int w = 0; w < current_subavalanches[q].size(); w++)
                            {
                                newsa.swaprecords_time_ID.push_back(pair<double, int>(swaprecords[nf][current_subavalanches[q][w]].time, current_subavalanches[q][w]));
                                swap_subavalancheID[nf][current_subavalanches[q][w]] = nsubavalanches;
                            }
                            subavalanches[nf].push_back(newsa);
                            nsubavalanches++;
                            current_subavalanches.erase(current_subavalanches.begin() + q);
                            isSubavStopped.erase(isSubavStopped.begin() + q);
                            tryErasing = true;
                            break;
                        }
                    }
                }
                //cout << " store+";
                //cout.flush();
            }
            // Store all lasting current subavalanches
            for (int q = 0; q < current_subavalanches.size(); q++)
            {
                subavalanche newsa;
                newsa.subavalancheID = nsubavalanches;
                for (int w = 0; w < current_subavalanches[q].size(); w++)
                {
                    newsa.swaprecords_time_ID.push_back(pair<double, int>(swaprecords[nf][current_subavalanches[q][w]].time, current_subavalanches[q][w]));
                    swap_subavalancheID[nf][current_subavalanches[q][w]] = nsubavalanches;
                }
                subavalanches[nf].push_back(newsa);
                nsubavalanches++;
            }

            // Calculate properties of subavalanches detected
            //cout << "\nsubavprop";
            //cout.flush();
            for (int nsa = 0; nsa < subavalanches[nf].size(); nsa++)
            {
                // Sort swaps by time
                sort(subavalanches[nf][nsa].swaprecords_time_ID.begin(), subavalanches[nf][nsa].swaprecords_time_ID.end());
                // Count
                subavalanches[nf][nsa].nswaps = subavalanches[nf][nsa].swaprecords_time_ID.size();
                // Timing
                subavalanches[nf][nsa].subavTimeStart = subavalanches[nf][nsa].swaprecords_time_ID[0].first;
                subavalanches[nf][nsa].subavTimeEnd = subavalanches[nf][nsa].swaprecords_time_ID[subavalanches[nf][nsa].nswaps - 1].first;
                subavalanches[nf][nsa].subavDuration = subavalanches[nf][nsa].subavTimeEnd - subavalanches[nf][nsa].subavTimeStart;
                int nswap_P50 = subavalanches[nf][nsa].nswaps / 2;
                if (subavalanches[nf][nsa].nswaps % 2 == 1) subavalanches[nf][nsa].subavTimeP50 = subavalanches[nf][nsa].swaprecords_time_ID[nswap_P50].first;
                else subavalanches[nf][nsa].subavTimeP50 = (subavalanches[nf][nsa].swaprecords_time_ID[nswap_P50].first + subavalanches[nf][nsa].swaprecords_time_ID[nswap_P50 - 1].first) * 0.50;
                // Geometry
                double d2max = 0.00;
                for (int ns1 = 0; ns1 < subavalanches[nf][nsa].nswaps; ns1++)
                {
                    int swapID1 = subavalanches[nf][nsa].swaprecords_time_ID[ns1].second;
                    for (int ns2 = ns1 + 1; ns2 < subavalanches[nf][nsa].nswaps; ns2++)
                    {
                        int swapID2 = subavalanches[nf][nsa].swaprecords_time_ID[ns2].second;
                        v3_math::spatial::vector<double> r = swaprecords[nf][swapID1].bondCenter - swaprecords[nf][swapID2].bondCenter;
                        double d2 = r.length2();
                        if (d2 > d2max) d2max = d2;
                    }
                }
                subavalanches[nf][nsa].subavDiameter = v3_math::msqrt(d2max) + avgBondLength;
                v3_math::spatial::vector<double> subav_center(0.00, 0.00, 0.00);
                subavalanches[nf][nsa].subavCenterX = subavalanches[nf][nsa].subavCenterY = 0.00;
                for (int ns1 = 0; ns1 < subavalanches[nf][nsa].nswaps; ns1++)
                {
                    int swapID1 = subavalanches[nf][nsa].swaprecords_time_ID[ns1].second;
                    subav_center += swaprecords[nf][swapID1].bondCenter;
                }
                subav_center *= 1.00 / subavalanches[nf][nsa].nswaps;
                subavalanches[nf][nsa].subavCenterX = subav_center.x;
                subavalanches[nf][nsa].subavCenterY = subav_center.y;
                subavalanches[nf][nsa].subavRg = 0.00;
                for (int ns1 = 0; ns1 < subavalanches[nf][nsa].nswaps; ns1++)
                {
                    int swapID1 = subavalanches[nf][nsa].swaprecords_time_ID[ns1].second;
                    v3_math::spatial::vector<double> dr = subav_center - swaprecords[nf][swapID1].bondCenter;
                    subavalanches[nf][nsa].subavRg += dr.length2();
                }
                subavalanches[nf][nsa].subavRg = v3_math::msqrt(subavalanches[nf][nsa].subavRg / double(subavalanches[nf][nsa].nswaps));
                if (subavalanches[nf][nsa].nswaps <= 1)
                    subavalanches[nf][nsa].subavI1 = subavalanches[nf][nsa].subavI2 = 0.00;
                else
                {
                    double Ixx = 0.00, Ixy = 0.00, Iyy = 0.00;
                    for (int ns1 = 0; ns1 < subavalanches[nf][nsa].nswaps; ns1++)
                    {
                        int swapID1 = subavalanches[nf][nsa].swaprecords_time_ID[ns1].second;
                        v3_math::spatial::vector<double> dr = subav_center - swaprecords[nf][swapID1].bondCenter;
                        Ixx += dr.y * dr.y;
                        Iyy += dr.x * dr.x;
                        Ixy -= dr.x * dr.y;
                    }
                    double d = v3_math::msqrt((Ixx - Iyy) * (Ixx - Iyy) + 4.00 * Ixy * Ixy);
                    subavalanches[nf][nsa].subavI1 = 0.5 * ((Ixx + Iyy) + d);
                    subavalanches[nf][nsa].subavI2 = 0.5 * ((Ixx + Iyy) - d);
                }
            }

            // Mark simultaneous subavalanches
            subavalanches_nsimultaneous[nf].resize(subavalanches[nf].size());
            for (int nsa = 0; nsa < subavalanches[nf].size(); nsa++)
                subavalanches_nsimultaneous[nf][nsa] = 0;
            for (int nsa = 0; nsa < subavalanches[nf].size(); nsa++)
            {
                for (int nsa2 = 0; nsa2 < subavalanches[nf].size(); nsa2++)
                {
                    if ((subavalanches[nf][nsa].subavTimeStart >= subavalanches[nf][nsa2].subavTimeStart &&
                        subavalanches[nf][nsa].subavTimeStart <= subavalanches[nf][nsa2].subavTimeEnd) ||
                        (subavalanches[nf][nsa].subavTimeEnd >= subavalanches[nf][nsa2].subavTimeStart &&
                            subavalanches[nf][nsa].subavTimeEnd <= subavalanches[nf][nsa2].subavTimeEnd) ||
                        (subavalanches[nf][nsa2].subavTimeStart >= subavalanches[nf][nsa].subavTimeStart &&
                            subavalanches[nf][nsa2].subavTimeStart <= subavalanches[nf][nsa].subavTimeEnd) ||
                        (subavalanches[nf][nsa2].subavTimeEnd >= subavalanches[nf][nsa].subavTimeStart &&
                            subavalanches[nf][nsa2].subavTimeEnd <= subavalanches[nf][nsa].subavTimeEnd))
                    {
                        subavalanches_nsimultaneous[nf][nsa]++;
                    }
                }
            }
        }
    } // if (subavalancheDetectionMechanism == 0)

    cout << "done.\n";
    /** / Detect subavalanches **/


/**
* This block combines subavalanches (cascades) into thermal avalanches
**/
    /** $ Detect avalanches **/
    cout << "Detect avalanches...\n";
    vector<vector<int>> swaprecords_avalancheID(nfiles);
    vector<vector<int>> subavalanche_avalancheID(nfiles);
    vector<vector<avalanche>> avalanches(nfiles);

    if (avalancheDetectionMechanism == 0)
    {
        cout << "avalancheDetectionMechanism == 0\n";

        for (int nf = 0; nf < nfiles; nf++)
        {
            swaprecords_avalancheID[nf].resize(swaprecords[nf].size());
            for (int q = 0; q < swaprecords[nf].size(); q++) swaprecords_avalancheID[nf][q] = -1;
            subavalanche_avalancheID[nf].resize(subavalanches[nf].size());
            for (int q = 0; q < subavalanches[nf].size(); q++) subavalanche_avalancheID[nf][q] = -1;

            // Sort subavalanches by timeP50
            vector<pair<double, int>> subavSorted_timeP50_subavID;
            for (int nsa = 0; nsa < subavalanches[nf].size(); nsa++)
                subavSorted_timeP50_subavID.push_back(pair<double, int> (subavalanches[nf][nsa].subavTimeP50, nsa));
            sort(subavSorted_timeP50_subavID.begin(), subavSorted_timeP50_subavID.end());

            // Detect interavalanche subavalanches
            vector<bool> subavSorted_timeP50_interavalanche(subavSorted_timeP50_subavID.size());
            if (subavSorted_timeP50_subavID.size() > 0) subavSorted_timeP50_interavalanche[0] = false;
            for (int q = 1; q < subavSorted_timeP50_subavID.size(); q++)
            {
                //int thisID = subavSorted_timeP50_subavID[q].second, prevID = subavSorted_timeP50_subavID[q - 1].second;
                double dttP50 = (subavSorted_timeP50_subavID[q].first - subavSorted_timeP50_subavID[q - 1].first) / subavSorted_timeP50_subavID[q].first;
                if (dttP50 >= avDetection_dtt_threshold) subavSorted_timeP50_interavalanche[q] = true;
                else subavSorted_timeP50_interavalanche[q] = false;
            }

            // Avalanche push_back cycle
            avalanche av;
            for (int q = 0; q < subavSorted_timeP50_subavID.size(); q++)
            {
                if (subavSorted_timeP50_interavalanche[q])
                {
                    // Store current avalanche 
                    av.avalancheID = avalanches[nf].size();
                    if (av.subavalanches_timeP50_ID.size() > 0)
                    {
                        sort(av.subavalanches_timeP50_ID.begin(), av.subavalanches_timeP50_ID.end());
                        sort(av.subavalanches_timeStart_ID.begin(), av.subavalanches_timeStart_ID.end());
                        sort(av.swaprecords_time_ID.begin(), av.swaprecords_time_ID.end());

                        for (int w = 0; w < av.swaprecords_time_ID.size(); w++)
                            swaprecords_avalancheID[nf][av.swaprecords_time_ID[w].second] = av.avalancheID;
                        for (int w = 0; w < av.subavalanches_timeP50_ID.size(); w++)
                            subavalanche_avalancheID[nf][av.subavalanches_timeP50_ID[w].second] = av.avalancheID;

                        av.avNsubavs = av.subavalanches_timeP50_ID.size();
                        av.avNswaps = av.swaprecords_time_ID.size();
                        av.avNswapsDir = 0;
                        for (int w = 0; w < av.swaprecords_time_ID.size(); w++)
                            av.avNswapsDir -= swaprecords[nf][av.swaprecords_time_ID[w].second].swapDirection;

                        av.avTimeStart = av.swaprecords_time_ID[0].first;
                        av.avTimeEnd = av.swaprecords_time_ID[av.swaprecords_time_ID.size() - 1].first;
                        av.avDuration = av.avTimeEnd - av.avTimeStart;
                        av.avThetaDuration = log10(av.avTimeEnd) - log10(av.avTimeStart);
                        av.avTimeAvg = 0.00;
                        av.avThetaAvg = 0.00;
                        v3_math::spatial::vector<double> avcenter(0.00, 0.00, 0.00);
                        for (int w = 0; w < av.swaprecords_time_ID.size(); w++)
                        {
                            av.avTimeAvg += av.swaprecords_time_ID[w].first;
                            av.avThetaAvg += log10(av.swaprecords_time_ID[w].first);
                            avcenter += swaprecords[nf][av.swaprecords_time_ID[w].second].bondCenter;
                        }
                        av.avTimeAvg /= double(av.swaprecords_time_ID.size());
                        av.avThetaAvg /= double(av.swaprecords_time_ID.size());
                        avcenter *= 1.00 / double(av.swaprecords_time_ID.size());
                        if (av.swaprecords_time_ID.size() % 2 == 1) av.avTimeP50 = av.swaprecords_time_ID[av.swaprecords_time_ID.size() / 2].first;
                        else av.avTimeP50 = (av.swaprecords_time_ID[av.swaprecords_time_ID.size() / 2].first + av.swaprecords_time_ID[av.swaprecords_time_ID.size() / 2 - 1].first) * 0.50;
                        av.avTg = 0.00;
                        av.avThetag = 0.00;
                        av.avRg = 0.00;
                        av.avDiameter = 0.00;
                        for (int w = 0; w < av.swaprecords_time_ID.size(); w++)
                        {
                            double dt = av.avTimeAvg - av.swaprecords_time_ID[w].first;
                            av.avTg += dt * dt;
                            double dtheta = av.avThetaAvg - log10(av.swaprecords_time_ID[w].first);
                            av.avThetag += dtheta * dtheta;
                            v3_math::spatial::vector<double> dr = avcenter - swaprecords[nf][av.swaprecords_time_ID[w].second].bondCenter;
                            av.avRg += dr.length2();
                            for (int e = w + 1; e < av.swaprecords_time_ID.size(); e++)
                            {
                                v3_math::spatial::vector<double> dr = swaprecords[nf][av.swaprecords_time_ID[w].second].bondCenter
                                    - swaprecords[nf][av.swaprecords_time_ID[e].second].bondCenter;
                                av.avDiameter = v3_math::max(av.avDiameter, dr.length());
                            }
                        }
                        av.avRg = v3_math::msqrt(av.avRg / double(av.avNswaps));
                        av.avTg = v3_math::msqrt(av.avTg / double(av.avNswaps));
                        av.avThetag = v3_math::msqrt(av.avThetag / double(av.avNswaps));

                        avalanches[nf].push_back(av);
                    }

                    // Start new avalanche
                    av.subavalanches_timeP50_ID.clear();
                    av.subavalanches_timeStart_ID.clear();
                    av.swaprecords_time_ID.clear();
                    av.subavalanches_timeP50_ID.push_back(subavSorted_timeP50_subavID[q]);
                    int subavID = subavSorted_timeP50_subavID[q].second;
                    av.subavalanches_timeStart_ID.push_back(pair<double, int>(subavalanches[nf][subavID].subavTimeStart, subavID));
                    for (int w = 0; w < subavalanches[nf][subavID].nswaps; w++)
                        av.swaprecords_time_ID.push_back(subavalanches[nf][subavID].swaprecords_time_ID[w]);
                }
                else
                {
                    // Add subav to the existing av
                    av.subavalanches_timeP50_ID.push_back(subavSorted_timeP50_subavID[q]);
                    int subavID = subavSorted_timeP50_subavID[q].second;
                    av.subavalanches_timeStart_ID.push_back(pair<double, int>(subavalanches[nf][subavID].subavTimeStart, subavID));
                    for (int w = 0; w < subavalanches[nf][subavID].nswaps; w++)
                        av.swaprecords_time_ID.push_back(subavalanches[nf][subavID].swaprecords_time_ID[w]);
                }
            }
            // Store the last avalanche 
            av.avalancheID = avalanches[nf].size();
            if (av.subavalanches_timeP50_ID.size() > 0)
            {
                sort(av.subavalanches_timeP50_ID.begin(), av.subavalanches_timeP50_ID.end());
                sort(av.subavalanches_timeStart_ID.begin(), av.subavalanches_timeStart_ID.end());
                sort(av.swaprecords_time_ID.begin(), av.swaprecords_time_ID.end());

                for (int w = 0; w < av.swaprecords_time_ID.size(); w++)
                    swaprecords_avalancheID[nf][av.swaprecords_time_ID[w].second] = av.avalancheID;
                for (int w = 0; w < av.subavalanches_timeP50_ID.size(); w++)
                    subavalanche_avalancheID[nf][av.subavalanches_timeP50_ID[w].second] = av.avalancheID;

                av.avNsubavs = av.subavalanches_timeP50_ID.size();
                av.avNswaps = av.swaprecords_time_ID.size();
                av.avNswapsDir = 0;
                for (int w = 0; w < av.swaprecords_time_ID.size(); w++)
                    av.avNswapsDir -= swaprecords[nf][av.swaprecords_time_ID[w].second].swapDirection;

                av.avTimeStart = av.swaprecords_time_ID[0].first;
                av.avTimeEnd = av.swaprecords_time_ID[av.swaprecords_time_ID.size() - 1].first;
                av.avDuration = av.avTimeEnd - av.avTimeStart;
                av.avThetaDuration = log10(av.avTimeEnd) - log10(av.avTimeStart);
                av.avTimeAvg = 0.00;
                av.avThetaAvg = 0.00;
                v3_math::spatial::vector<double> avcenter(0.00, 0.00, 0.00);
                for (int w = 0; w < av.swaprecords_time_ID.size(); w++)
                {
                    av.avTimeAvg += av.swaprecords_time_ID[w].first;
                    av.avThetaAvg += log10(av.swaprecords_time_ID[w].first);
                    avcenter += swaprecords[nf][av.swaprecords_time_ID[w].second].bondCenter;
                }
                av.avTimeAvg /= double(av.swaprecords_time_ID.size());
                av.avThetaAvg /= double(av.swaprecords_time_ID.size());
                avcenter *= 1.00 / double(av.swaprecords_time_ID.size());
                if (av.swaprecords_time_ID.size() % 2 == 1) av.avTimeP50 = av.swaprecords_time_ID[av.swaprecords_time_ID.size() / 2].first;
                else av.avTimeP50 = (av.swaprecords_time_ID[av.swaprecords_time_ID.size() / 2].first + av.swaprecords_time_ID[av.swaprecords_time_ID.size() / 2 - 1].first) * 0.50;
                av.avTg = 0.00;
                av.avThetag = 0.00;
                av.avRg = 0.00;
                av.avDiameter = 0.00;
                for (int w = 0; w < av.swaprecords_time_ID.size(); w++)
                {
                    double dt = av.avTimeAvg - av.swaprecords_time_ID[w].first;
                    av.avTg += dt * dt;
                    double dtheta = av.avThetaAvg - log10(av.swaprecords_time_ID[w].first);
                    av.avThetag += dtheta * dtheta;
                    v3_math::spatial::vector<double> dr = avcenter - swaprecords[nf][av.swaprecords_time_ID[w].second].bondCenter;
                    av.avRg += dr.length2();
                    for (int e = w + 1; e < av.swaprecords_time_ID.size(); e++)
                    {
                        v3_math::spatial::vector<double> dr = swaprecords[nf][av.swaprecords_time_ID[w].second].bondCenter
                            - swaprecords[nf][av.swaprecords_time_ID[e].second].bondCenter;
                        av.avDiameter = v3_math::max(av.avDiameter, dr.length());
                    }
                }
                av.avRg = v3_math::msqrt(av.avRg / double(av.avNswaps));
                av.avTg = v3_math::msqrt(av.avTg / double(av.avNswaps));
                av.avThetag = v3_math::msqrt(av.avThetag / double(av.avNswaps));

                avalanches[nf].push_back(av);
            }
        } // for (int nf = 0; nf < nfiles; nf++)
    } // if (avalancheDetectionMechanism == 0)
    else if (avalancheDetectionMechanism == 1)
    {
        cout << "avalancheDetectionMechanism == 1\n";

        for (int nf = 0; nf < nfiles; nf++)
        {
            swaprecords_avalancheID[nf].resize(swaprecords[nf].size());
            for (int q = 0; q < swaprecords[nf].size(); q++) swaprecords_avalancheID[nf][q] = -1;
            subavalanche_avalancheID[nf].resize(subavalanches[nf].size());
            for (int q = 0; q < subavalanches[nf].size(); q++) subavalanche_avalancheID[nf][q] = -1;

            // Sort subavalanches by timeP50
            vector<pair<double, int>> subavSorted_timeP50_subavID;
            for (int nsa = 0; nsa < subavalanches[nf].size(); nsa++)
                subavSorted_timeP50_subavID.push_back(pair<double, int>(subavalanches[nf][nsa].subavTimeP50, nsa));
            sort(subavSorted_timeP50_subavID.begin(), subavSorted_timeP50_subavID.end());

            // Local groups
            int navalanches = 0;
            vector<double> subav_avalanche_dr(subavalanches[nf].size());
            vector<vector<int>> current_avalanches;
            for (int ns = 0; ns < subavalanches[nf].size(); ns++)
            {
                // Calculate dtt and dr for each of the current avalanches
                int subavID = subavSorted_timeP50_subavID[ns].second;
                int nav = current_avalanches.size();
                vector<double> dtP50_last(nav), dttP50_last(nav), dr_threshold_min(nav), dttP50_threshold_min(nav);
                vector<bool> threshold_min_defined(nav);
                for (int q = 0; q < nav; q++)
                {
                    int avSize = current_avalanches[q].size();
                    int lastSubavID = current_avalanches[q][avSize - 1];
                    dtP50_last[q] = subavSorted_timeP50_subavID[ns].first - subavalanches[nf][lastSubavID].subavTimeP50;
                    dttP50_last[q] = dtP50_last[q] / subavSorted_timeP50_subavID[ns].first;
                    threshold_min_defined[q] = false;
                    for (int w = 0; w < avSize; w++)
                    {
                        int wID = current_avalanches[q][avSize - 1 - w];
                        double dtP50 = subavSorted_timeP50_subavID[ns].first - subavalanches[nf][wID].subavTimeP50;
                        double dttP50 = dtP50 / subavSorted_timeP50_subavID[ns].first;
                        if (dttP50 > avDetection_dtt_threshold) break;
                        double dr2 = (subavalanches[nf][subavID].subavCenterX - subavalanches[nf][wID].subavCenterX) * (subavalanches[nf][subavID].subavCenterX - subavalanches[nf][wID].subavCenterX) +
                            (subavalanches[nf][subavID].subavCenterY - subavalanches[nf][wID].subavCenterY) * (subavalanches[nf][subavID].subavCenterY - subavalanches[nf][wID].subavCenterY);
                        double dr = v3_math::msqrt(dr2);
                        if (threshold_min_defined[q])
                        {
                            if (dr_threshold_min[q] > dr)
                            {
                                dr_threshold_min[q] = dr;
                                dttP50_threshold_min[q] = dttP50;
                            }
                        }
                        else
                        {
                            dr_threshold_min[q] = dr;
                            dttP50_threshold_min[q] = dttP50;
                            threshold_min_defined[q] = true;
                        }
                    }
                }

                // Determine the closest avalanche
                int closestAv = -1;
                double dr_min = 0.00;
                for (int q = 0; q < nav; q++)
                {
                    if (threshold_min_defined[q])
                    {
                        if (dr_threshold_min[q] < avDetection_dr_threshold)
                        {
                            if (closestAv < 0)
                            {
                                dr_min = dr_threshold_min[q];
                                closestAv = q;
                            }
                            else if (dr_min > dr_threshold_min[q])
                            {
                                dr_min = dr_threshold_min[q];
                                closestAv = q;
                            }
                        }
                    }
                }

                // Determine the distance-to-previous associated with avalanche
                subav_avalanche_dr[ns] = 0.00;
                if (ns > 0)
                {
                    int prevSubavID = subavSorted_timeP50_subavID[ns - 1].second;
                    double dr2 = (subavalanches[nf][subavID].subavCenterX - subavalanches[nf][prevSubavID].subavCenterX) * (subavalanches[nf][subavID].subavCenterX - subavalanches[nf][prevSubavID].subavCenterX) +
                        (subavalanches[nf][subavID].subavCenterY - subavalanches[nf][prevSubavID].subavCenterY) * (subavalanches[nf][subavID].subavCenterY - subavalanches[nf][prevSubavID].subavCenterY);
                    double dr = v3_math::msqrt(dr2);
                    subav_avalanche_dr[ns] = dr;
                }
                if (closestAv >= 0) subav_avalanche_dr[ns] = dr_min;

                // Check if the new subavalanche breaks any of the avalanche sequences
                vector<bool> isAvStopped(nav);
                int nstop = 0;
                for (int q = 0; q < nav; q++)
                {
                    isAvStopped[q] = false;
                    if (dttP50_last[q] > avDetection_dtt_threshold)
                    {
                        isAvStopped[q] = true;
                        nstop++;
                        continue;
                    }
                }

                // Add the subavalanche to the closest avalanche (if exists) or to the new avalanche
                if (closestAv >= 0)
                {
                    current_avalanches[closestAv].push_back(subavSorted_timeP50_subavID[ns].second);
                }
                else
                {
                    vector<int> newa(1);
                    newa[0] = subavSorted_timeP50_subavID[ns].second;
                    current_avalanches.push_back(newa);
                    isAvStopped.push_back(false);
                }

                // Store & forget stopped current avalanches
                bool tryErasing = true;
                while (tryErasing)
                {
                    tryErasing = false;
                    for (int q = 0; q < current_avalanches.size(); q++)
                    {
                        if (isAvStopped[q])
                        {
                            //cout << "\nStore subav " << q << " of " << current_subavalanches.size()
                            //    << " with " << current_subavalanches[q].size() << " events";
                            //cout.flush();
                            avalanche newav;
                            newav.avalancheID = navalanches;
                            for (int w = 0; w < current_avalanches[q].size(); w++)
                            {
                                int subavID = current_avalanches[q][w];
                                newav.subavalanches_timeP50_ID.push_back(pair<double,int> (subavalanches[nf][subavID].subavTimeP50, subavID));
                                newav.subavalanches_timeStart_ID.push_back(pair<double, int>(subavalanches[nf][subavID].subavTimeStart, subavID));
                                subavalanche_avalancheID[nf][subavID] = newav.avalancheID;
                                for (int w = 0; w < subavalanches[nf][subavID].nswaps; w++)
                                {
                                    newav.swaprecords_time_ID.push_back(subavalanches[nf][subavID].swaprecords_time_ID[w]);
                                    swaprecords_avalancheID[nf][subavalanches[nf][subavID].swaprecords_time_ID[w].second] = newav.avalancheID;
                                }
                            }
                            avalanches[nf].push_back(newav);
                            navalanches++;
                            current_avalanches.erase(current_avalanches.begin() + q);
                            isAvStopped.erase(isAvStopped.begin() + q);
                            tryErasing = true;
                            break;
                        }
                    }
                }
            }

            // Store all lasting current subavalanches
            for (int q = 0; q < current_avalanches.size(); q++)
            {
                avalanche newav;
                newav.avalancheID = navalanches;
                for (int w = 0; w < current_avalanches[q].size(); w++)
                {
                    int subavID = current_avalanches[q][w];
                    newav.subavalanches_timeP50_ID.push_back(pair<double, int>(subavalanches[nf][subavID].subavTimeP50, subavID));
                    newav.subavalanches_timeStart_ID.push_back(pair<double, int>(subavalanches[nf][subavID].subavTimeStart, subavID));
                    subavalanche_avalancheID[nf][subavID] = newav.avalancheID;
                    for (int w = 0; w < subavalanches[nf][subavID].nswaps; w++)
                    {
                        newav.swaprecords_time_ID.push_back(subavalanches[nf][subavID].swaprecords_time_ID[w]);
                        swaprecords_avalancheID[nf][subavalanches[nf][subavID].swaprecords_time_ID[w].second] = newav.avalancheID;
                    }
                }
                avalanches[nf].push_back(newav);
                navalanches++;
            }

            // Calculate avalanche properties
            for (int nav = 0; nav < avalanches[nf].size(); nav++)
            {
                avalanches[nf][nav].avNsubavs = avalanches[nf][nav].subavalanches_timeP50_ID.size();
                avalanches[nf][nav].avNswaps = avalanches[nf][nav].swaprecords_time_ID.size();
                avalanches[nf][nav].avNswapsDir = 0;
                for (int w = 0; w < avalanches[nf][nav].swaprecords_time_ID.size(); w++)
                    avalanches[nf][nav].avNswapsDir -= swaprecords[nf][avalanches[nf][nav].swaprecords_time_ID[w].second].swapDirection;

                avalanches[nf][nav].avTimeStart = avalanches[nf][nav].swaprecords_time_ID[0].first;
                avalanches[nf][nav].avTimeEnd = avalanches[nf][nav].swaprecords_time_ID[avalanches[nf][nav].swaprecords_time_ID.size() - 1].first;
                avalanches[nf][nav].avDuration = avalanches[nf][nav].avTimeEnd - avalanches[nf][nav].avTimeStart;
                avalanches[nf][nav].avThetaDuration = log10(avalanches[nf][nav].avTimeEnd) - log10(avalanches[nf][nav].avTimeStart);
                avalanches[nf][nav].avTimeAvg = 0.00;
                avalanches[nf][nav].avThetaAvg = 0.00;
                v3_math::spatial::vector<double> avcenter(0.00, 0.00, 0.00);
                for (int w = 0; w < avalanches[nf][nav].swaprecords_time_ID.size(); w++)
                {
                    avalanches[nf][nav].avTimeAvg += avalanches[nf][nav].swaprecords_time_ID[w].first;
                    avalanches[nf][nav].avThetaAvg += log10(avalanches[nf][nav].swaprecords_time_ID[w].first);
                    avcenter += swaprecords[nf][avalanches[nf][nav].swaprecords_time_ID[w].second].bondCenter;
                }
                avalanches[nf][nav].avTimeAvg /= double(avalanches[nf][nav].swaprecords_time_ID.size());
                avalanches[nf][nav].avThetaAvg /= double(avalanches[nf][nav].swaprecords_time_ID.size());
                avcenter *= 1.00 / double(avalanches[nf][nav].swaprecords_time_ID.size());
                if (avalanches[nf][nav].swaprecords_time_ID.size() % 2 == 1) avalanches[nf][nav].avTimeP50 = avalanches[nf][nav].swaprecords_time_ID[avalanches[nf][nav].swaprecords_time_ID.size() / 2].first;
                else avalanches[nf][nav].avTimeP50 = (avalanches[nf][nav].swaprecords_time_ID[avalanches[nf][nav].swaprecords_time_ID.size() / 2].first + avalanches[nf][nav].swaprecords_time_ID[avalanches[nf][nav].swaprecords_time_ID.size() / 2 - 1].first) * 0.50;
                avalanches[nf][nav].avTg = 0.00;
                avalanches[nf][nav].avThetag = 0.00;
                avalanches[nf][nav].avRg = 0.00;
                avalanches[nf][nav].avDiameter = 0.00;
                for (int w = 0; w < avalanches[nf][nav].swaprecords_time_ID.size(); w++)
                {
                    double dt = avalanches[nf][nav].avTimeAvg - avalanches[nf][nav].swaprecords_time_ID[w].first;
                    avalanches[nf][nav].avTg += dt * dt;
                    double dtheta = avalanches[nf][nav].avThetaAvg - log10(avalanches[nf][nav].swaprecords_time_ID[w].first);
                    avalanches[nf][nav].avThetag += dtheta * dtheta;
                    v3_math::spatial::vector<double> dr = avcenter - swaprecords[nf][avalanches[nf][nav].swaprecords_time_ID[w].second].bondCenter;
                    avalanches[nf][nav].avRg += dr.length2();
                    for (int e = w + 1; e < avalanches[nf][nav].swaprecords_time_ID.size(); e++)
                    {
                        v3_math::spatial::vector<double> dr = swaprecords[nf][avalanches[nf][nav].swaprecords_time_ID[w].second].bondCenter
                            - swaprecords[nf][avalanches[nf][nav].swaprecords_time_ID[e].second].bondCenter;
                        avalanches[nf][nav].avDiameter = v3_math::max(avalanches[nf][nav].avDiameter, dr.length());
                    }
                }
                avalanches[nf][nav].avRg = v3_math::msqrt(avalanches[nf][nav].avRg / double(avalanches[nf][nav].avNswaps));
                avalanches[nf][nav].avTg = v3_math::msqrt(avalanches[nf][nav].avTg / double(avalanches[nf][nav].avNswaps));
                avalanches[nf][nav].avThetag = v3_math::msqrt(avalanches[nf][nav].avThetag / double(avalanches[nf][nav].avNswaps));
            }
        }
    }
    else if (avalancheDetectionMechanism == 2)
    {
        cout << "avalancheDetectionMechanism == 2\n";

        for (int nf = 0; nf < nfiles; nf++)
        {
            swaprecords_avalancheID[nf].resize(swaprecords[nf].size());
            for (int q = 0; q < swaprecords[nf].size(); q++) swaprecords_avalancheID[nf][q] = -1;
            subavalanche_avalancheID[nf].resize(subavalanches[nf].size());
            for (int q = 0; q < subavalanches[nf].size(); q++) subavalanche_avalancheID[nf][q] = -1;

            // Sort subavalanches by timeP50
            vector<pair<double, int>> subavSorted_timeP50_subavID;
            for (int nsa = 0; nsa < subavalanches[nf].size(); nsa++)
                subavSorted_timeP50_subavID.push_back(pair<double, int>(subavalanches[nf][nsa].subavTimeP50, nsa));
            sort(subavSorted_timeP50_subavID.begin(), subavSorted_timeP50_subavID.end());

            // Local groups
            int navalanches = 0;
            vector<vector<int>> current_avalanches_subavID;
            vector<int> current_avalanches_nswaps;

            // Scan all subavalanches subsequently
            for (int ns = 0; ns < subavalanches[nf].size(); ns++)
            {
                // Check what avalanches it can fit into
                int scanSubavID = subavSorted_timeP50_subavID[ns].second;
                int ncav = current_avalanches_subavID.size();
                vector<bool> canjoin(ncav);
                for (int na = 0; na < ncav; na++)
                {
                    canjoin[na] = false;
                    int nsubav = current_avalanches_subavID[na].size();
                    for (int q = nsubav - 1; (q >= 0 && !canjoin[na]); q--)
                    {
                        int prevSubavID = current_avalanches_subavID[na][q];
                        if (subavalanches[nf][scanSubavID].nswaps < avDetection_subavNswaps_threshold) continue;
                        double dtP50 = subavalanches[nf][scanSubavID].subavTimeP50 - subavalanches[nf][prevSubavID].subavTimeP50;
                        double dttP50 = dtP50 / subavalanches[nf][scanSubavID].subavTimeP50;
                        if (dttP50 > avDetection_dtt_threshold) break;
                        double dr2 = (subavalanches[nf][scanSubavID].subavCenterX - subavalanches[nf][prevSubavID].subavCenterX) * (subavalanches[nf][scanSubavID].subavCenterX - subavalanches[nf][prevSubavID].subavCenterX) +
                            (subavalanches[nf][scanSubavID].subavCenterY - subavalanches[nf][prevSubavID].subavCenterY) * (subavalanches[nf][scanSubavID].subavCenterY - subavalanches[nf][prevSubavID].subavCenterY);
                        double dr = v3_math::msqrt(dr2);
                        if (dr > avDetection_dr_threshold) continue;
                        canjoin[na] = true;
                    }
                }

                // Choose the biggest avalanche to fit
                int biggestAvNum = -1, biggestAvSize = -1;
                for (int na = 0; na < ncav; na++)
                {
                    if (canjoin[na])
                    {
                        if (biggestAvSize < current_avalanches_nswaps[na])
                        {
                            biggestAvSize = current_avalanches_nswaps[na];
                            biggestAvNum = na;
                        }
                    }
                }

                // Add the subavalanche to the closest avalanche (if exists) or to the new avalanche
                if (biggestAvNum >= 0)
                {
                    current_avalanches_subavID[biggestAvNum].push_back(scanSubavID);
                    current_avalanches_nswaps[biggestAvNum] += subavalanches[nf][scanSubavID].nswaps;
                }
                else
                {
                    vector<int> newa(1);
                    newa[0] = scanSubavID;
                    current_avalanches_subavID.push_back(newa);
                    current_avalanches_nswaps.push_back(subavalanches[nf][scanSubavID].nswaps);
                }
            }

            // Store all lasting current subavalanches
            for (int q = 0; q < current_avalanches_subavID.size(); q++)
            {
                avalanche newav;
                newav.avalancheID = navalanches;
                for (int w = 0; w < current_avalanches_subavID[q].size(); w++)
                {
                    int subavID = current_avalanches_subavID[q][w];
                    newav.subavalanches_timeP50_ID.push_back(pair<double, int>(subavalanches[nf][subavID].subavTimeP50, subavID));
                    newav.subavalanches_timeStart_ID.push_back(pair<double, int>(subavalanches[nf][subavID].subavTimeStart, subavID));
                    subavalanche_avalancheID[nf][subavID] = newav.avalancheID;
                    for (int w = 0; w < subavalanches[nf][subavID].nswaps; w++)
                    {
                        newav.swaprecords_time_ID.push_back(subavalanches[nf][subavID].swaprecords_time_ID[w]);
                        swaprecords_avalancheID[nf][subavalanches[nf][subavID].swaprecords_time_ID[w].second] = newav.avalancheID;
                    }
                }
                avalanches[nf].push_back(newav);
                navalanches++;
            }

            // Calculate avalanche properties
            for (int nav = 0; nav < avalanches[nf].size(); nav++)
            {
                avalanches[nf][nav].avNsubavs = avalanches[nf][nav].subavalanches_timeP50_ID.size();
                avalanches[nf][nav].avNswaps = avalanches[nf][nav].swaprecords_time_ID.size();
                avalanches[nf][nav].avNswapsDir = 0;
                for (int w = 0; w < avalanches[nf][nav].swaprecords_time_ID.size(); w++)
                    avalanches[nf][nav].avNswapsDir -= swaprecords[nf][avalanches[nf][nav].swaprecords_time_ID[w].second].swapDirection;

                avalanches[nf][nav].avTimeStart = avalanches[nf][nav].swaprecords_time_ID[0].first;
                avalanches[nf][nav].avTimeEnd = avalanches[nf][nav].swaprecords_time_ID[avalanches[nf][nav].swaprecords_time_ID.size() - 1].first;
                avalanches[nf][nav].avDuration = avalanches[nf][nav].avTimeEnd - avalanches[nf][nav].avTimeStart;
                avalanches[nf][nav].avThetaDuration = log10(avalanches[nf][nav].avTimeEnd) - log10(avalanches[nf][nav].avTimeStart);
                avalanches[nf][nav].avTimeAvg = 0.00;
                avalanches[nf][nav].avThetaAvg = 0.00;
                v3_math::spatial::vector<double> avcenter(0.00, 0.00, 0.00);
                for (int w = 0; w < avalanches[nf][nav].swaprecords_time_ID.size(); w++)
                {
                    avalanches[nf][nav].avTimeAvg += avalanches[nf][nav].swaprecords_time_ID[w].first;
                    avalanches[nf][nav].avThetaAvg += log10(avalanches[nf][nav].swaprecords_time_ID[w].first);
                    avcenter += swaprecords[nf][avalanches[nf][nav].swaprecords_time_ID[w].second].bondCenter;
                }
                avalanches[nf][nav].avTimeAvg /= double(avalanches[nf][nav].swaprecords_time_ID.size());
                avalanches[nf][nav].avThetaAvg /= double(avalanches[nf][nav].swaprecords_time_ID.size());
                avcenter *= 1.00 / double(avalanches[nf][nav].swaprecords_time_ID.size());
                if (avalanches[nf][nav].swaprecords_time_ID.size() % 2 == 1) avalanches[nf][nav].avTimeP50 = avalanches[nf][nav].swaprecords_time_ID[avalanches[nf][nav].swaprecords_time_ID.size() / 2].first;
                else avalanches[nf][nav].avTimeP50 = (avalanches[nf][nav].swaprecords_time_ID[avalanches[nf][nav].swaprecords_time_ID.size() / 2].first + avalanches[nf][nav].swaprecords_time_ID[avalanches[nf][nav].swaprecords_time_ID.size() / 2 - 1].first) * 0.50;
                avalanches[nf][nav].avTg = 0.00;
                avalanches[nf][nav].avThetag = 0.00;
                avalanches[nf][nav].avRg = 0.00;
                avalanches[nf][nav].avDiameter = 0.00;
                for (int w = 0; w < avalanches[nf][nav].swaprecords_time_ID.size(); w++)
                {
                    double dt = avalanches[nf][nav].avTimeAvg - avalanches[nf][nav].swaprecords_time_ID[w].first;
                    avalanches[nf][nav].avTg += dt * dt;
                    double dtheta = avalanches[nf][nav].avThetaAvg - log10(avalanches[nf][nav].swaprecords_time_ID[w].first);
                    avalanches[nf][nav].avThetag += dtheta * dtheta;
                    v3_math::spatial::vector<double> dr = avcenter - swaprecords[nf][avalanches[nf][nav].swaprecords_time_ID[w].second].bondCenter;
                    avalanches[nf][nav].avRg += dr.length2();
                    for (int e = w + 1; e < avalanches[nf][nav].swaprecords_time_ID.size(); e++)
                    {
                        v3_math::spatial::vector<double> dr = swaprecords[nf][avalanches[nf][nav].swaprecords_time_ID[w].second].bondCenter
                            - swaprecords[nf][avalanches[nf][nav].swaprecords_time_ID[e].second].bondCenter;
                        avalanches[nf][nav].avDiameter = v3_math::max(avalanches[nf][nav].avDiameter, dr.length());
                    }
                }
                avalanches[nf][nav].avRg = v3_math::msqrt(avalanches[nf][nav].avRg / double(avalanches[nf][nav].avNswaps));
                avalanches[nf][nav].avTg = v3_math::msqrt(avalanches[nf][nav].avTg / double(avalanches[nf][nav].avNswaps));
                avalanches[nf][nav].avThetag = v3_math::msqrt(avalanches[nf][nav].avThetag / double(avalanches[nf][nav].avNswaps));
            }
        }
    }
    else if (avalancheDetectionMechanism == 3)
    {
        cout << "avalancheDetectionMechanism == 3\n";

        for (int nf = 0; nf < nfiles; nf++)
        {
            swaprecords_avalancheID[nf].resize(swaprecords[nf].size());
            for (int q = 0; q < swaprecords[nf].size(); q++) swaprecords_avalancheID[nf][q] = -1;
            subavalanche_avalancheID[nf].resize(subavalanches[nf].size());
            for (int q = 0; q < subavalanches[nf].size(); q++) subavalanche_avalancheID[nf][q] = -1;

            // Sort subavalanches by timeP50
            vector<pair<double, int>> subavSorted_timeP50_subavID;
            for (int nsa = 0; nsa < subavalanches[nf].size(); nsa++)
                subavSorted_timeP50_subavID.push_back(pair<double, int>(subavalanches[nf][nsa].subavTimeP50, nsa));
            sort(subavSorted_timeP50_subavID.begin(), subavSorted_timeP50_subavID.end());

            // Local groups
            int navalanches = 0;
            vector<vector<int>> current_avalanches_subavID;
            vector<int> current_avalanches_nswaps;

            // Scan all subavalanches subsequently
            for (int ns = 0; ns < subavalanches[nf].size(); ns++)
            {
                // Check what avalanches it can fit into
                int scanSubavID = subavSorted_timeP50_subavID[ns].second;
                int ncav = current_avalanches_subavID.size();
                vector<bool> canjoin(ncav);
                for (int na = 0; na < ncav; na++)
                {
                    canjoin[na] = false;
                    int nsubav = current_avalanches_subavID[na].size();
                    for (int q = nsubav - 1; (q >= 0 && !canjoin[na]); q--)
                    {
                        int prevSubavID = current_avalanches_subavID[na][q];
                        if (subavalanches[nf][scanSubavID].nswaps < avDetection_subavNswaps_threshold) continue;
                        double dtP50 = subavalanches[nf][scanSubavID].subavTimeP50 - subavalanches[nf][prevSubavID].subavTimeP50;
                        double dttP50 = dtP50 / subavalanches[nf][scanSubavID].subavTimeP50;
                        if (dtP50 > avDetection_dt_minimum && dttP50 > avDetection_dtt_threshold) break;
                        double dr2 = (subavalanches[nf][scanSubavID].subavCenterX - subavalanches[nf][prevSubavID].subavCenterX) * (subavalanches[nf][scanSubavID].subavCenterX - subavalanches[nf][prevSubavID].subavCenterX) +
                            (subavalanches[nf][scanSubavID].subavCenterY - subavalanches[nf][prevSubavID].subavCenterY) * (subavalanches[nf][scanSubavID].subavCenterY - subavalanches[nf][prevSubavID].subavCenterY);
                        double dr = v3_math::msqrt(dr2);
                        if (dr > avDetection_dr_threshold) continue;
                        canjoin[na] = true;
                    }
                }

                // Choose the biggest avalanche to fit
                int biggestAvNum = -1, biggestAvSize = -1;
                for (int na = 0; na < ncav; na++)
                {
                    if (canjoin[na])
                    {
                        if (biggestAvSize < current_avalanches_nswaps[na])
                        {
                            biggestAvSize = current_avalanches_nswaps[na];
                            biggestAvNum = na;
                        }
                    }
                }

                // Add the subavalanche to the closest avalanche (if exists) or to the new avalanche
                if (biggestAvNum >= 0)
                {
                    current_avalanches_subavID[biggestAvNum].push_back(scanSubavID);
                    current_avalanches_nswaps[biggestAvNum] += subavalanches[nf][scanSubavID].nswaps;
                }
                else
                {
                    vector<int> newa(1);
                    newa[0] = scanSubavID;
                    current_avalanches_subavID.push_back(newa);
                    current_avalanches_nswaps.push_back(subavalanches[nf][scanSubavID].nswaps);
                }
            }

            // Store all lasting current subavalanches
            for (int q = 0; q < current_avalanches_subavID.size(); q++)
            {
                avalanche newav;
                newav.avalancheID = navalanches;
                for (int w = 0; w < current_avalanches_subavID[q].size(); w++)
                {
                    int subavID = current_avalanches_subavID[q][w];
                    newav.subavalanches_timeP50_ID.push_back(pair<double, int>(subavalanches[nf][subavID].subavTimeP50, subavID));
                    newav.subavalanches_timeStart_ID.push_back(pair<double, int>(subavalanches[nf][subavID].subavTimeStart, subavID));
                    subavalanche_avalancheID[nf][subavID] = newav.avalancheID;
                    for (int w = 0; w < subavalanches[nf][subavID].nswaps; w++)
                    {
                        newav.swaprecords_time_ID.push_back(subavalanches[nf][subavID].swaprecords_time_ID[w]);
                        swaprecords_avalancheID[nf][subavalanches[nf][subavID].swaprecords_time_ID[w].second] = newav.avalancheID;
                    }
                }
                avalanches[nf].push_back(newav);
                navalanches++;
            }

            // Calculate avalanche properties
            for (int nav = 0; nav < avalanches[nf].size(); nav++)
            {
                avalanches[nf][nav].avNsubavs = avalanches[nf][nav].subavalanches_timeP50_ID.size();
                avalanches[nf][nav].avNswaps = avalanches[nf][nav].swaprecords_time_ID.size();
                avalanches[nf][nav].avNswapsDir = 0;
                for (int w = 0; w < avalanches[nf][nav].swaprecords_time_ID.size(); w++)
                    avalanches[nf][nav].avNswapsDir -= swaprecords[nf][avalanches[nf][nav].swaprecords_time_ID[w].second].swapDirection;

                avalanches[nf][nav].avTimeStart = avalanches[nf][nav].swaprecords_time_ID[0].first;
                avalanches[nf][nav].avTimeEnd = avalanches[nf][nav].swaprecords_time_ID[avalanches[nf][nav].swaprecords_time_ID.size() - 1].first;
                avalanches[nf][nav].avDuration = avalanches[nf][nav].avTimeEnd - avalanches[nf][nav].avTimeStart;
                avalanches[nf][nav].avThetaDuration = log10(avalanches[nf][nav].avTimeEnd) - log10(avalanches[nf][nav].avTimeStart);
                avalanches[nf][nav].avTimeAvg = 0.00;
                avalanches[nf][nav].avThetaAvg = 0.00;
                v3_math::spatial::vector<double> avcenter(0.00, 0.00, 0.00);
                for (int w = 0; w < avalanches[nf][nav].swaprecords_time_ID.size(); w++)
                {
                    avalanches[nf][nav].avTimeAvg += avalanches[nf][nav].swaprecords_time_ID[w].first;
                    avalanches[nf][nav].avThetaAvg += log10(avalanches[nf][nav].swaprecords_time_ID[w].first);
                    avcenter += swaprecords[nf][avalanches[nf][nav].swaprecords_time_ID[w].second].bondCenter;
                }
                avalanches[nf][nav].avTimeAvg /= double(avalanches[nf][nav].swaprecords_time_ID.size());
                avalanches[nf][nav].avThetaAvg /= double(avalanches[nf][nav].swaprecords_time_ID.size());
                avcenter *= 1.00 / double(avalanches[nf][nav].swaprecords_time_ID.size());
                if (avalanches[nf][nav].swaprecords_time_ID.size() % 2 == 1) avalanches[nf][nav].avTimeP50 = avalanches[nf][nav].swaprecords_time_ID[avalanches[nf][nav].swaprecords_time_ID.size() / 2].first;
                else avalanches[nf][nav].avTimeP50 = (avalanches[nf][nav].swaprecords_time_ID[avalanches[nf][nav].swaprecords_time_ID.size() / 2].first + avalanches[nf][nav].swaprecords_time_ID[avalanches[nf][nav].swaprecords_time_ID.size() / 2 - 1].first) * 0.50;
                avalanches[nf][nav].avTg = 0.00;
                avalanches[nf][nav].avThetag = 0.00;
                avalanches[nf][nav].avRg = 0.00;
                avalanches[nf][nav].avDiameter = 0.00;
                for (int w = 0; w < avalanches[nf][nav].swaprecords_time_ID.size(); w++)
                {
                    double dt = avalanches[nf][nav].avTimeAvg - avalanches[nf][nav].swaprecords_time_ID[w].first;
                    avalanches[nf][nav].avTg += dt * dt;
                    double dtheta = avalanches[nf][nav].avThetaAvg - log10(avalanches[nf][nav].swaprecords_time_ID[w].first);
                    avalanches[nf][nav].avThetag += dtheta * dtheta;
                    v3_math::spatial::vector<double> dr = avcenter - swaprecords[nf][avalanches[nf][nav].swaprecords_time_ID[w].second].bondCenter;
                    avalanches[nf][nav].avRg += dr.length2();
                    for (int e = w + 1; e < avalanches[nf][nav].swaprecords_time_ID.size(); e++)
                    {
                        v3_math::spatial::vector<double> dr = swaprecords[nf][avalanches[nf][nav].swaprecords_time_ID[w].second].bondCenter
                            - swaprecords[nf][avalanches[nf][nav].swaprecords_time_ID[e].second].bondCenter;
                        avalanches[nf][nav].avDiameter = v3_math::max(avalanches[nf][nav].avDiameter, dr.length());
                    }
                }
                avalanches[nf][nav].avRg = v3_math::msqrt(avalanches[nf][nav].avRg / double(avalanches[nf][nav].avNswaps));
                avalanches[nf][nav].avTg = v3_math::msqrt(avalanches[nf][nav].avTg / double(avalanches[nf][nav].avNswaps));
                avalanches[nf][nav].avThetag = v3_math::msqrt(avalanches[nf][nav].avThetag / double(avalanches[nf][nav].avNswaps));
            }
        }
    }
    else if (avalancheDetectionMechanism == 4)
    {
        cout << "avalancheDetectionMechanism == 4\n";

        for (int nf = 0; nf < nfiles; nf++)
        {
            swaprecords_avalancheID[nf].resize(swaprecords[nf].size());
            for (int q = 0; q < swaprecords[nf].size(); q++) swaprecords_avalancheID[nf][q] = -1;
            subavalanche_avalancheID[nf].resize(subavalanches[nf].size());
            for (int q = 0; q < subavalanches[nf].size(); q++) subavalanche_avalancheID[nf][q] = -1;

            // Sort subavalanches by timeP50
            vector<pair<double, int>> subavSorted_timeP50_subavID;
            for (int nsa = 0; nsa < subavalanches[nf].size(); nsa++)
                subavSorted_timeP50_subavID.push_back(pair<double, int>(subavalanches[nf][nsa].subavTimeP50, nsa));
            sort(subavSorted_timeP50_subavID.begin(), subavSorted_timeP50_subavID.end());

            // Local groups
            int navalanches = 0;
            vector<vector<int>> current_avalanches_subavID;
            vector<int> current_avalanches_nswaps;

            // Scan all subavalanches subsequently
            for (int ns = 0; ns < subavalanches[nf].size(); ns++)
            {
                // Check what avalanches it can fit into
                int scanSubavID = subavSorted_timeP50_subavID[ns].second;
                int ncav = current_avalanches_subavID.size();
                vector<bool> canjoin(ncav);
                for (int na = 0; na < ncav; na++)
                {
                    canjoin[na] = false;
                    int nsubav = current_avalanches_subavID[na].size();
                    for (int q = nsubav - 1; (q >= 0 && !canjoin[na]); q--)
                    {
                        int prevSubavID = current_avalanches_subavID[na][q];
                        if (subavalanches[nf][scanSubavID].nswaps < avDetection_subavNswaps_threshold) continue;
                        double dtP50 = subavalanches[nf][scanSubavID].subavTimeP50 - subavalanches[nf][prevSubavID].subavTimeP50;
                        double dttP50 = dtP50 / subavalanches[nf][scanSubavID].subavTimeP50;
                        if (dtP50 > avDetection_dt_minimum && dttP50 > avDetection_dtt_threshold) break;
                        double dr2 = (subavalanches[nf][scanSubavID].subavCenterX - subavalanches[nf][prevSubavID].subavCenterX) * (subavalanches[nf][scanSubavID].subavCenterX - subavalanches[nf][prevSubavID].subavCenterX) +
                            (subavalanches[nf][scanSubavID].subavCenterY - subavalanches[nf][prevSubavID].subavCenterY) * (subavalanches[nf][scanSubavID].subavCenterY - subavalanches[nf][prevSubavID].subavCenterY);
                        double dr = v3_math::msqrt(dr2);
                        if (dr > dtP50 * avDetection_velocity_threshold) continue;
                        canjoin[na] = true;
                    }
                }

                // Choose the biggest avalanche to fit
                int biggestAvNum = -1, biggestAvSize = -1;
                for (int na = 0; na < ncav; na++)
                {
                    if (canjoin[na])
                    {
                        if (biggestAvSize < current_avalanches_nswaps[na])
                        {
                            biggestAvSize = current_avalanches_nswaps[na];
                            biggestAvNum = na;
                        }
                    }
                }

                // Add the subavalanche to the closest avalanche (if exists) or to the new avalanche
                if (biggestAvNum >= 0)
                {
                    current_avalanches_subavID[biggestAvNum].push_back(scanSubavID);
                    current_avalanches_nswaps[biggestAvNum] += subavalanches[nf][scanSubavID].nswaps;
                }
                else
                {
                    vector<int> newa(1);
                    newa[0] = scanSubavID;
                    current_avalanches_subavID.push_back(newa);
                    current_avalanches_nswaps.push_back(subavalanches[nf][scanSubavID].nswaps);
                }
            }

            // Store all lasting current subavalanches
            for (int q = 0; q < current_avalanches_subavID.size(); q++)
            {
                avalanche newav;
                newav.avalancheID = navalanches;
                for (int w = 0; w < current_avalanches_subavID[q].size(); w++)
                {
                    int subavID = current_avalanches_subavID[q][w];
                    newav.subavalanches_timeP50_ID.push_back(pair<double, int>(subavalanches[nf][subavID].subavTimeP50, subavID));
                    newav.subavalanches_timeStart_ID.push_back(pair<double, int>(subavalanches[nf][subavID].subavTimeStart, subavID));
                    subavalanche_avalancheID[nf][subavID] = newav.avalancheID;
                    for (int w = 0; w < subavalanches[nf][subavID].nswaps; w++)
                    {
                        newav.swaprecords_time_ID.push_back(subavalanches[nf][subavID].swaprecords_time_ID[w]);
                        swaprecords_avalancheID[nf][subavalanches[nf][subavID].swaprecords_time_ID[w].second] = newav.avalancheID;
                    }
                }
                avalanches[nf].push_back(newav);
                navalanches++;
            }

            // Calculate avalanche properties
            for (int nav = 0; nav < avalanches[nf].size(); nav++)
            {
                avalanches[nf][nav].avNsubavs = avalanches[nf][nav].subavalanches_timeP50_ID.size();
                avalanches[nf][nav].avNswaps = avalanches[nf][nav].swaprecords_time_ID.size();
                avalanches[nf][nav].avNswapsDir = 0;
                for (int w = 0; w < avalanches[nf][nav].swaprecords_time_ID.size(); w++)
                    avalanches[nf][nav].avNswapsDir -= swaprecords[nf][avalanches[nf][nav].swaprecords_time_ID[w].second].swapDirection;

                avalanches[nf][nav].avTimeStart = avalanches[nf][nav].swaprecords_time_ID[0].first;
                avalanches[nf][nav].avTimeEnd = avalanches[nf][nav].swaprecords_time_ID[avalanches[nf][nav].swaprecords_time_ID.size() - 1].first;
                avalanches[nf][nav].avDuration = avalanches[nf][nav].avTimeEnd - avalanches[nf][nav].avTimeStart;
                avalanches[nf][nav].avTimeAvg = 0.00;
                v3_math::spatial::vector<double> avcenter(0.00, 0.00, 0.00);
                for (int w = 0; w < avalanches[nf][nav].swaprecords_time_ID.size(); w++)
                {
                    avalanches[nf][nav].avTimeAvg += avalanches[nf][nav].swaprecords_time_ID[w].first;
                    avcenter += swaprecords[nf][avalanches[nf][nav].swaprecords_time_ID[w].second].bondCenter;
                }
                avalanches[nf][nav].avTimeAvg /= double(avalanches[nf][nav].swaprecords_time_ID.size());
                avcenter *= 1.00 / double(avalanches[nf][nav].swaprecords_time_ID.size());
                if (avalanches[nf][nav].swaprecords_time_ID.size() % 2 == 1) avalanches[nf][nav].avTimeP50 = avalanches[nf][nav].swaprecords_time_ID[avalanches[nf][nav].swaprecords_time_ID.size() / 2].first;
                else avalanches[nf][nav].avTimeP50 = (avalanches[nf][nav].swaprecords_time_ID[avalanches[nf][nav].swaprecords_time_ID.size() / 2].first + avalanches[nf][nav].swaprecords_time_ID[avalanches[nf][nav].swaprecords_time_ID.size() / 2 - 1].first) * 0.50;
                avalanches[nf][nav].avTg = 0.00;
                avalanches[nf][nav].avRg = 0.00;
                avalanches[nf][nav].avDiameter = 0.00;
                for (int w = 0; w < avalanches[nf][nav].swaprecords_time_ID.size(); w++)
                {
                    double dt = avalanches[nf][nav].avTimeAvg - avalanches[nf][nav].swaprecords_time_ID[w].first;
                    avalanches[nf][nav].avTg += dt * dt;
                    v3_math::spatial::vector<double> dr = avcenter - swaprecords[nf][avalanches[nf][nav].swaprecords_time_ID[w].second].bondCenter;
                    avalanches[nf][nav].avRg += dr.length2();
                    for (int e = w + 1; e < avalanches[nf][nav].swaprecords_time_ID.size(); e++)
                    {
                        v3_math::spatial::vector<double> dr = swaprecords[nf][avalanches[nf][nav].swaprecords_time_ID[w].second].bondCenter
                            - swaprecords[nf][avalanches[nf][nav].swaprecords_time_ID[e].second].bondCenter;
                        avalanches[nf][nav].avDiameter = v3_math::max(avalanches[nf][nav].avDiameter, dr.length());
                    }
                }
                avalanches[nf][nav].avRg = v3_math::msqrt(avalanches[nf][nav].avRg / double(avalanches[nf][nav].avNswaps));
                avalanches[nf][nav].avTg = v3_math::msqrt(avalanches[nf][nav].avTg / double(avalanches[nf][nav].avNswaps));
            }
        }
    }
    else if (avalancheDetectionMechanism == 5)
    {
        cout << "avalancheDetectionMechanism == 5\n";

        for (int nf = 0; nf < nfiles; nf++)
        {
            swaprecords_avalancheID[nf].resize(swaprecords[nf].size());
            for (int q = 0; q < swaprecords[nf].size(); q++) swaprecords_avalancheID[nf][q] = -1;
            subavalanche_avalancheID[nf].resize(subavalanches[nf].size());
            for (int q = 0; q < subavalanches[nf].size(); q++) subavalanche_avalancheID[nf][q] = -1;

            // Sort subavalanches by timeP50
            vector<pair<double, int>> subavSorted_timeP50_subavID;
            for (int nsa = 0; nsa < subavalanches[nf].size(); nsa++)
                subavSorted_timeP50_subavID.push_back(pair<double, int>(subavalanches[nf][nsa].subavTimeP50, nsa));
            sort(subavSorted_timeP50_subavID.begin(), subavSorted_timeP50_subavID.end());

            // Local groups
            int navalanches = 0;
            vector<pair<int, vector<int>>> current_avalanches_nswaps_subavID;
            vector<vector<int>> subavID_possibleAvalanches(subavalanches[nf].size());

            // Scan all subavalanches subsequently
            for (int ns = 0; ns < subavSorted_timeP50_subavID.size(); ns++)
            {
                // Check what avalanches it can fit into
                int scanSubavID = subavSorted_timeP50_subavID[ns].second;
                int ncav = current_avalanches_nswaps_subavID.size();
                bool nojoin = true;
                vector<bool> canjoin(ncav);
                for (int na = 0; na < ncav; na++)
                {
                    canjoin[na] = false;
                    int nsubav = current_avalanches_nswaps_subavID[na].second.size();
                    for (int q = nsubav - 1; (q >= 0 && !canjoin[na]); q--)
                    {
                        int prevSubavID = current_avalanches_nswaps_subavID[na].second[q];
                        double dtP50 = subavalanches[nf][scanSubavID].subavTimeP50 - subavalanches[nf][prevSubavID].subavTimeP50;
                        double dttP50 = dtP50 / subavalanches[nf][scanSubavID].subavTimeP50;
                        if (dttP50 > avDetection_dtt_threshold) break;
                        double dr2 = (subavalanches[nf][scanSubavID].subavCenterX - subavalanches[nf][prevSubavID].subavCenterX) * (subavalanches[nf][scanSubavID].subavCenterX - subavalanches[nf][prevSubavID].subavCenterX) +
                            (subavalanches[nf][scanSubavID].subavCenterY - subavalanches[nf][prevSubavID].subavCenterY) * (subavalanches[nf][scanSubavID].subavCenterY - subavalanches[nf][prevSubavID].subavCenterY);
                        double dr = v3_math::msqrt(dr2);
                        if (dr > dtP50 * avDetection_velocity_threshold) continue;
                        canjoin[na] = true;
                        nojoin = false;
                    }
                }

                // Add the subavalanche all possible avalanches, or form new
                if (nojoin) 
                {
                    pair<int, vector<int>> newa;
                    newa.first = subavalanches[nf][scanSubavID].nswaps;
                    newa.second.push_back(scanSubavID);
                    current_avalanches_nswaps_subavID.push_back(newa);
                    subavID_possibleAvalanches[scanSubavID].push_back(current_avalanches_nswaps_subavID.size() - 1);
                }
                else
                {
                    for (int na = 0; na < ncav; na++)
                    {
                        if (canjoin[na])
                        {
                            current_avalanches_nswaps_subavID[na].first += subavalanches[nf][scanSubavID].nswaps;
                            current_avalanches_nswaps_subavID[na].second.push_back(scanSubavID);
                            subavID_possibleAvalanches[scanSubavID].push_back(na);
                        }
                    }
                }
            }

            // Decide what goes where, the biggest possible avalanche wins it all
            vector<vector<int>> biggest_avalanches_subavID;
            while (current_avalanches_nswaps_subavID.size() > 0)
            {
                // Sort so the biggest avalanche is the last one
                sort(current_avalanches_nswaps_subavID.begin(), current_avalanches_nswaps_subavID.end());

                // Push the biggest avalanche
                int nbiggest = current_avalanches_nswaps_subavID.size() - 1;
                biggest_avalanches_subavID.push_back(current_avalanches_nswaps_subavID[nbiggest].second);

                // Remove subavalanches of the biggest avalanche from other avalanches
                for (int nsa = 0; nsa < current_avalanches_nswaps_subavID[nbiggest].second.size(); nsa++)
                {
                    int subavID = current_avalanches_nswaps_subavID[nbiggest].second[nsa];
                    for (int q = 0; q < current_avalanches_nswaps_subavID.size(); q++)
                    {
                        vector<int>::iterator it = std::find(current_avalanches_nswaps_subavID[q].second.begin(), current_avalanches_nswaps_subavID[q].second.end(), subavID);
                        if (it != current_avalanches_nswaps_subavID[q].second.end())
                        { // Found => erase
                            current_avalanches_nswaps_subavID[q].second.erase(it);
                            current_avalanches_nswaps_subavID[q].first -= subavalanches[nf][subavID].nswaps;
                        }
                    }
                }

                // Check the biggest avalanche is now zero
                if (current_avalanches_nswaps_subavID[nbiggest].first == 0 && current_avalanches_nswaps_subavID[nbiggest].second.size() == 0)
                    cout << "+";
                else 
                    cout << "!";

                // Erase the biggest avalanche from the vector
                current_avalanches_nswaps_subavID.erase(current_avalanches_nswaps_subavID.begin() + nbiggest);
            }
            
            // Store biggest subavalanches
            for (int q = 0; q < biggest_avalanches_subavID.size(); q++)
            {
                avalanche newav;
                newav.avalancheID = navalanches;
                for (int w = 0; w < biggest_avalanches_subavID[q].size(); w++)
                {
                    int subavID = biggest_avalanches_subavID[q][w];
                    newav.subavalanches_timeP50_ID.push_back(pair<double, int>(subavalanches[nf][subavID].subavTimeP50, subavID));
                    newav.subavalanches_timeStart_ID.push_back(pair<double, int>(subavalanches[nf][subavID].subavTimeStart, subavID));
                    subavalanche_avalancheID[nf][subavID] = newav.avalancheID;
                    for (int w = 0; w < subavalanches[nf][subavID].nswaps; w++)
                    {
                        newav.swaprecords_time_ID.push_back(subavalanches[nf][subavID].swaprecords_time_ID[w]);
                        swaprecords_avalancheID[nf][subavalanches[nf][subavID].swaprecords_time_ID[w].second] = newav.avalancheID;
                    }
                }
                avalanches[nf].push_back(newav);
                navalanches++;
            }

            // Calculate avalanche properties
            for (int nav = 0; nav < avalanches[nf].size(); nav++)
            {
                avalanches[nf][nav].avNsubavs = avalanches[nf][nav].subavalanches_timeP50_ID.size();
                avalanches[nf][nav].avNswaps = avalanches[nf][nav].swaprecords_time_ID.size();
                avalanches[nf][nav].avNswapsDir = 0;
                for (int w = 0; w < avalanches[nf][nav].swaprecords_time_ID.size(); w++)
                    avalanches[nf][nav].avNswapsDir -= swaprecords[nf][avalanches[nf][nav].swaprecords_time_ID[w].second].swapDirection;

                avalanches[nf][nav].avTimeStart = avalanches[nf][nav].swaprecords_time_ID[0].first;
                avalanches[nf][nav].avTimeEnd = avalanches[nf][nav].swaprecords_time_ID[avalanches[nf][nav].swaprecords_time_ID.size() - 1].first;
                avalanches[nf][nav].avDuration = avalanches[nf][nav].avTimeEnd - avalanches[nf][nav].avTimeStart;
                avalanches[nf][nav].avTimeAvg = 0.00;
                v3_math::spatial::vector<double> avcenter(0.00, 0.00, 0.00);
                for (int w = 0; w < avalanches[nf][nav].swaprecords_time_ID.size(); w++)
                {
                    avalanches[nf][nav].avTimeAvg += avalanches[nf][nav].swaprecords_time_ID[w].first;
                    avcenter += swaprecords[nf][avalanches[nf][nav].swaprecords_time_ID[w].second].bondCenter;
                }
                avalanches[nf][nav].avTimeAvg /= double(avalanches[nf][nav].swaprecords_time_ID.size());
                avcenter *= 1.00 / double(avalanches[nf][nav].swaprecords_time_ID.size());
                if (avalanches[nf][nav].swaprecords_time_ID.size() % 2 == 1) avalanches[nf][nav].avTimeP50 = avalanches[nf][nav].swaprecords_time_ID[avalanches[nf][nav].swaprecords_time_ID.size() / 2].first;
                else avalanches[nf][nav].avTimeP50 = (avalanches[nf][nav].swaprecords_time_ID[avalanches[nf][nav].swaprecords_time_ID.size() / 2].first + avalanches[nf][nav].swaprecords_time_ID[avalanches[nf][nav].swaprecords_time_ID.size() / 2 - 1].first) * 0.50;
                avalanches[nf][nav].avTg = 0.00;
                avalanches[nf][nav].avRg = 0.00;
                avalanches[nf][nav].avDiameter = 0.00;
                for (int w = 0; w < avalanches[nf][nav].swaprecords_time_ID.size(); w++)
                {
                    double dt = avalanches[nf][nav].avTimeAvg - avalanches[nf][nav].swaprecords_time_ID[w].first;
                    avalanches[nf][nav].avTg += dt * dt;
                    v3_math::spatial::vector<double> dr = avcenter - swaprecords[nf][avalanches[nf][nav].swaprecords_time_ID[w].second].bondCenter;
                    avalanches[nf][nav].avRg += dr.length2();
                    for (int e = w + 1; e < avalanches[nf][nav].swaprecords_time_ID.size(); e++)
                    {
                        v3_math::spatial::vector<double> dr = swaprecords[nf][avalanches[nf][nav].swaprecords_time_ID[w].second].bondCenter
                            - swaprecords[nf][avalanches[nf][nav].swaprecords_time_ID[e].second].bondCenter;
                        avalanches[nf][nav].avDiameter = v3_math::max(avalanches[nf][nav].avDiameter, dr.length());
                    }
                }
                avalanches[nf][nav].avRg = v3_math::msqrt(avalanches[nf][nav].avRg / double(avalanches[nf][nav].avNswaps));
                avalanches[nf][nav].avTg = v3_math::msqrt(avalanches[nf][nav].avTg / double(avalanches[nf][nav].avNswaps));
            }
        }
    }

    cout << "done.\n";
    /** / Detect avalanches **/


/**
* This block calculates statistics on bond snaps
**/
    /** $ Collect swap statistics **/
    cout << "\nCollect swap statistics...\n";
    // ALL SWAPS
    cout << "All swaps\n";
    vector<v3_math::planar::rectangle<double>> boxLimits(nfiles);
    v3_math::Histogramm<double> swap_epoch(0.00, swapStat_nepoch, swapStat_nepoch* swapStat_nepochfrac, v3_math::HISTMODE_EXT);
    vector<v3_math::Histogramm<double>> swap_logdt_epoch(swapStat_nepoch* swapStat_nepochfrac + 1);
    vector<v3_math::AvgErr<double>> swap_avgdt_epoch(swapStat_nepoch* swapStat_nepochfrac + 1);
    vector<v3_math::Histogramm<double>> swap_dtheta_epoch(swapStat_nepoch* swapStat_nepochfrac + 1);
    vector<v3_math::CHistogramm<double>> swap_cdtheta_epoch(swapStat_nepoch* swapStat_nepochfrac + 1);
    vector<v3_math::CHistogramm<double>> swap_cdthetaav_epoch(swapStat_nepoch* swapStat_nepochfrac + 1);
    vector<v3_math::CHistogramm<double>> swap_cdtt_epoch(swapStat_nepoch* swapStat_nepochfrac + 1);
    vector<v3_math::CHistogramm<double>> swap_cdttav_epoch(swapStat_nepoch* swapStat_nepochfrac + 1);
    v3_math::Histogramm<double> histdr(swapStat_dr_min, swapStat_dr_max, swapStat_dr_nbins, v3_math::HISTMODE_CUTOFF);
    v3_math::Histogramm<double> histlogdr(swapStat_logdr_min, swapStat_logdr_max, swapStat_logdr_nbins, v3_math::HISTMODE_CUTOFF);
    v3_math::Histogramm2D<double> histdtdr(swapStat_dtdr_logdt_min, swapStat_dtdr_logdt_max, swapStat_dtdr_logdr_min,
        swapStat_dtdr_logdr_max, swapStat_dtdr_logdt_nbins, swapStat_dtdr_logdr_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm2D<double> histdtmindr(swapStat_dtdr_logdt_min, swapStat_dtdr_logdt_max, swapStat_dtdr_logdr_min,
        swapStat_dtdr_logdr_max, swapStat_dtdr_logdt_nbins, swapStat_dtdr_logdr_nbins, v3_math::HISTMODE_EXT);
    vector<double> dthetaBinBorders;
    dthetaBinBorders.push_back(0.00);
    dthetaBinBorders.push_back(0.000000001);
    dthetaBinBorders.push_back(0.0000000015);
    dthetaBinBorders.push_back(0.000000002);
    dthetaBinBorders.push_back(0.000000003);
    dthetaBinBorders.push_back(0.0000000045);
    dthetaBinBorders.push_back(0.0000000065);
    dthetaBinBorders.push_back(0.00000001);
    dthetaBinBorders.push_back(0.000000015);
    dthetaBinBorders.push_back(0.00000002);
    dthetaBinBorders.push_back(0.00000003);
    dthetaBinBorders.push_back(0.000000045);
    dthetaBinBorders.push_back(0.000000065);
    dthetaBinBorders.push_back(0.0000001);
    dthetaBinBorders.push_back(0.00000015);
    dthetaBinBorders.push_back(0.0000002);
    dthetaBinBorders.push_back(0.0000003);
    dthetaBinBorders.push_back(0.00000045);
    dthetaBinBorders.push_back(0.00000065);
    dthetaBinBorders.push_back(0.000001);
    dthetaBinBorders.push_back(0.0000015);
    dthetaBinBorders.push_back(0.000002);
    dthetaBinBorders.push_back(0.000003);
    dthetaBinBorders.push_back(0.0000045);
    dthetaBinBorders.push_back(0.0000065);
    dthetaBinBorders.push_back(0.00001);
    dthetaBinBorders.push_back(0.000015);
    dthetaBinBorders.push_back(0.00002);
    dthetaBinBorders.push_back(0.00003);
    dthetaBinBorders.push_back(0.000045);
    dthetaBinBorders.push_back(0.000065);
    dthetaBinBorders.push_back(0.0001);
    dthetaBinBorders.push_back(0.00015);
    dthetaBinBorders.push_back(0.0002);
    dthetaBinBorders.push_back(0.0003);
    dthetaBinBorders.push_back(0.00045);
    dthetaBinBorders.push_back(0.00065);
    dthetaBinBorders.push_back(0.001);
    dthetaBinBorders.push_back(0.0015);
    dthetaBinBorders.push_back(0.002);
    dthetaBinBorders.push_back(0.003);
    dthetaBinBorders.push_back(0.0045);
    dthetaBinBorders.push_back(0.0065);
    dthetaBinBorders.push_back(0.01);
    dthetaBinBorders.push_back(0.015);
    dthetaBinBorders.push_back(0.02);
    dthetaBinBorders.push_back(0.03);
    dthetaBinBorders.push_back(0.045);
    dthetaBinBorders.push_back(0.065);
    dthetaBinBorders.push_back(0.1);
    dthetaBinBorders.push_back(0.2);
    dthetaBinBorders.push_back(0.3);
    dthetaBinBorders.push_back(0.45);
    dthetaBinBorders.push_back(0.65);
    dthetaBinBorders.push_back(1.0);
    dthetaBinBorders.push_back(2.0);
    v3_math::Histogramm2D<double> hist_ij_logdtdr(swapStat_ij_dtdr_logdt_min, swapStat_ij_dtdr_logdt_max, 
        0.00, pow(10.00, swapStat_ij_dtdr_logdr_max), 
        swapStat_ij_dtdr_logdt_nbins, swapStat_ij_dtdr_logdr_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm2D<double> hist_ij_logdtlogdr(swapStat_ij_dtdr_logdt_min, swapStat_ij_dtdr_logdt_max, swapStat_ij_dtdr_logdr_min,
        swapStat_ij_dtdr_logdr_max, swapStat_ij_dtdr_logdt_nbins, swapStat_ij_dtdr_logdr_nbins, v3_math::HISTMODE_EXT);
    v3_math::WHistogramm<double> swapStat_ijCount_whistdt(swapStat_ij_dt_min, swapStat_ij_dt_max, swapStat_ij_dt_nbins, v3_math::HISTMODE_EXT),
        swapStat_ijNoweight_whistdt(swapStat_ij_dt_min, swapStat_ij_dt_max, swapStat_ij_dt_nbins, v3_math::HISTMODE_EXT);
    v3_math::WHistogramm<double> swapStat_ijCount_whistdr(swapStat_ij_dr_min, swapStat_ij_dr_max, swapStat_ij_dr_nbins, v3_math::HISTMODE_EXT),
        swapStat_ijNoweight_whistdr(swapStat_ij_dr_min, swapStat_ij_dr_max, swapStat_ij_dr_nbins, v3_math::HISTMODE_EXT);
    v3_math::WHistogramm<double> swapStat_ijCount_whistdr_alltime(swapStat_ij_dr_min, swapStat_ij_dr_max, swapStat_ij_dr_nbins, v3_math::HISTMODE_EXT),
        swapStat_ijNoweight_whistdr_alltime(swapStat_ij_dr_min, swapStat_ij_dr_max, swapStat_ij_dr_nbins, v3_math::HISTMODE_EXT);
    vector<v3_math::WHistogramm<double>> swapStat_ijCount_whistdtdr(swapStat_ij_dr_nbins), swapStat_ijNoweight_whistdtdr(swapStat_ij_dr_nbins);
    double swapStat_ijCount_rbinsize = (swapStat_ij_dr_max - swapStat_ij_dr_min) / double(swapStat_ij_dr_nbins);
    for (int q = 0; q < swapStat_ij_dr_nbins; q++)
    {
        swapStat_ijCount_whistdtdr[q].change(swapStat_ij_dt_min, swapStat_ij_dt_max, swapStat_ij_dt_nbins, v3_math::HISTMODE_CUTOFF);
        swapStat_ijNoweight_whistdtdr[q].change(swapStat_ij_dt_min, swapStat_ij_dt_max, swapStat_ij_dt_nbins, v3_math::HISTMODE_CUTOFF);
    }
    v3_math::WHistogramm<double> swapStat_ijCount_whistdtheta(swapStat_ij_dtheta_min, swapStat_ij_dtheta_max, swapStat_ij_dtheta_nbins, v3_math::HISTMODE_EXT),
        swapStat_ijNoweight_whistdtheta(swapStat_ij_dtheta_min, swapStat_ij_dtheta_max, swapStat_ij_dtheta_nbins, v3_math::HISTMODE_EXT);
    vector<v3_math::WHistogramm<double>> swapStat_ijCount_whistdthetadr(swapStat_ij_dr_nbins), swapStat_ijNoweight_whistdthetadr(swapStat_ij_dr_nbins);
    for (int q = 0; q < swapStat_ij_dr_nbins; q++)
    {
        swapStat_ijCount_whistdthetadr[q].change(swapStat_ij_dtheta_min, swapStat_ij_dtheta_max, swapStat_ij_dtheta_nbins, v3_math::HISTMODE_CUTOFF);
        swapStat_ijNoweight_whistdthetadr[q].change(swapStat_ij_dtheta_min, swapStat_ij_dtheta_max, swapStat_ij_dtheta_nbins, v3_math::HISTMODE_CUTOFF);
    }
    for (int q = 0; q < swapStat_nepoch * swapStat_nepochfrac + 1; q++)
    {
        swap_logdt_epoch[q].change(-3.00, 7.00, 10 * 3, v3_math::HISTMODE_EXT);
        swap_dtheta_epoch[q].change(0.00, 0.50, 100, v3_math::HISTMODE_EXT);
        swap_cdtheta_epoch[q].change(dthetaBinBorders, v3_math::CHISTMODE_EXT);
        swap_cdthetaav_epoch[q].change(dthetaBinBorders, v3_math::CHISTMODE_EXT);
        swap_cdtt_epoch[q].change(dthetaBinBorders, v3_math::CHISTMODE_EXT);
        swap_cdttav_epoch[q].change(dthetaBinBorders, v3_math::CHISTMODE_EXT);
    }
    vector<int> nswaps_logt;
    nswaps_logt.push_back(100);
    nswaps_logt.push_back(200);
    nswaps_logt.push_back(500);
    nswaps_logt.push_back(1000);
    nswaps_logt.push_back(2000);
    nswaps_logt.push_back(3000);
    nswaps_logt.push_back(4000);
    nswaps_logt.push_back(5000);
    nswaps_logt.push_back(10000);
    vector<v3_math::Histogramm<double>> swapStat_nswaps_logthist(nswaps_logt.size());
    for (int q = 0; q < nswaps_logt.size(); q++) swapStat_nswaps_logthist[q].change(0.00, 8.00, 8 * 2);
    for (int nf = 0; nf < nfiles; nf++)
    {
        cout << "\nFile " << nf + 1 << " / " << nfiles;
        cout.flush();
        for (int ns = 1; ns < swaprecords[nf].size(); ns++)
        {
            for (int q = 0; q < nswaps_logt.size(); q++)
            {
                if (ns == nswaps_logt[q] && swaprecords[nf][ns].time > 0.00)
                {
                    swapStat_nswaps_logthist[q].add(log10(swaprecords[nf][ns].time));
                    if (log10(swaprecords[nf][ns].time) > 6.50)
                        cout << "log10(t): nswaps=" << ns << ", t=" << swaprecords[nf][ns].time << ", log10(t)=" << log10(swaprecords[nf][ns].time) << "\n";
                }
            }

            if (swaprecords[nf][ns].time < analysis_time_min) continue;
            if (analysis_time_max > 0.00 && swaprecords[nf][ns].time > analysis_time_max) continue;

            double swapTime_log10 = log10(swaprecords[nf][ns].time);
            int ne = swapTime_log10 - swapStat_epochStart_log10;
            int nef = (swapTime_log10 - swapStat_epochStart_log10) * double(swapStat_nepochfrac);
            swap_epoch.add(swapTime_log10 - swapStat_epochStart_log10);
            double swap_dt = swaprecords[nf][ns].time - swaprecords[nf][ns - 1].time;
            if (swap_dt > 0.00)
            {
                double swap_dtt = swap_dt / swaprecords[nf][ns].time;
                double log10dt = log10(swap_dt);
                double swap_dr = (swaprecords[nf][ns].bondCenter - swaprecords[nf][ns - 1].bondCenter).length();
                swap_logdt_epoch[0].add(log10dt);
                swap_avgdt_epoch[0].add(swap_dt);
                swap_dtheta_epoch[0].add(swap_dtt);
                swap_cdtheta_epoch[0].add(swap_dtt);
                swap_cdtt_epoch[0].add(swap_dtt);
                if (nef >= 0 && nef < swapStat_nepoch * swapStat_nepochfrac)
                {
                    swap_logdt_epoch[nef + 1].add(log10dt);
                    swap_avgdt_epoch[nef + 1].add(swap_dt);
                    swap_dtheta_epoch[nef + 1].add(swap_dtt);
                    swap_cdtheta_epoch[nef + 1].add(swap_dtt);
                    swap_cdtt_epoch[nef + 1].add(swap_dtt);
                }
                histdr.add(swap_dr);
                if (swap_dr > 0.00)
                {
                    histdtdr.add(log10dt, log10(swap_dr));
                    histlogdr.add(log10(swap_dr));
                }
            }
        } // for (int ns = 1; ns < swaprecords[nf].size(); ns++)

        boxLimits[nf].L.x = system_posLeft_x;
        boxLimits[nf].L.y = system_posLeft_y;
        boxLimits[nf].R.x = system_posLeft_x + swaprecords[nf][0].xsize;
        boxLimits[nf].R.y = system_posLeft_y + swaprecords[nf][0].ysize;
        cout << "box limits: x " << boxLimits[nf].L.x << " to " << boxLimits[nf].R.x << ", y " << boxLimits[nf].L.y << " to " << boxLimits[nf].R.y << "\n";
        for (int nsi = 0; nsi < swaprecords[nf].size(); nsi++)
        {
            if (swaprecords[nf][nsi].time < analysis_time_min) continue;
            if (analysis_time_max > 0.00 && swaprecords[nf][nsi].time > analysis_time_max) continue;

            double tmax = v3_math::max(analysis_time_max, maxsimtime[nf]);
            double dtmax = tmax - swaprecords[nf][nsi].time;
            double dthetamax = log10(tmax) - log10(swaprecords[nf][nsi].time);
            swapStat_ijNoweight_whistdt.addAllBinsBefore(log10(dtmax), 1.00);
            swapStat_ijNoweight_whistdtheta.addAllBinsBefore(log10(dthetamax), 1.00);
            for (int rbin = 0; rbin < swapStat_ij_dr_nbins; rbin++)
            {
                swapStat_ijNoweight_whistdtdr[rbin].addAllBinsBefore(log10(dtmax), 1.00);
                swapStat_ijNoweight_whistdthetadr[rbin].addAllBinsBefore(log10(dthetamax), 1.00);
            }
            v3_math::planar::point<double> boxLeft(system_posLeft_x, system_posLeft_y), 
                boxRight(system_posLeft_x + swaprecords[nf][nsi].xsize, system_posLeft_y + swaprecords[nf][nsi].ysize);
            v3_math::planar::rectangle<double> box(boxLeft, boxRight);
            bool update = false;
            if (boxLimits[nf].L.x > box.L.x) {boxLimits[nf].L.x = box.L.x; update = true;}
            if (boxLimits[nf].L.y > box.L.y) {boxLimits[nf].L.y = box.L.y; update = true;}
            if (boxLimits[nf].R.x < box.R.x) {boxLimits[nf].R.x = box.R.x; update = true;}
            if (boxLimits[nf].R.y < box.R.y) {boxLimits[nf].R.y = box.R.y; update = true;}
            if (update)
                cout << "nsi = " << nsi << ", new box limits: x " << boxLimits[nf].L.x << " to " << boxLimits[nf].R.x << ", y " << boxLimits[nf].L.y << " to " << boxLimits[nf].R.y << "\n";
            v3_math::planar::point<double> bondCenter(swaprecords[nf][nsi].bondCenter.x, swaprecords[nf][nsi].bondCenter.y);
            for (int rbin = 0; rbin < swapStat_ij_dr_nbins; rbin++)
            {
                double rin = swapStat_ij_dr_min + double(rbin) * swapStat_ijCount_rbinsize;
                double rout = swapStat_ij_dr_min + double(rbin + 1) * swapStat_ijCount_rbinsize;
                v3_math::planar::circle<double> Cin(bondCenter, rin), Cout(bondCenter, rout);
                double Sin = v3_math::planar::intersectarea(box, Cin);
                double Sout = v3_math::planar::intersectarea(box, Cout);
                double Sbin = Sout - Sin;
                swapStat_ijNoweight_whistdr.add((rin+rout)*0.50, Sbin);
                swapStat_ijNoweight_whistdr_alltime.add((rin+rout)*0.50, Sbin);
            }

            for (int nsj = nsi + 1; nsj < swaprecords[nf].size(); nsj++)
            {
                if (swaprecords[nf][nsj].time < analysis_time_min) continue;
                if (analysis_time_max > 0.00 && swaprecords[nf][nsj].time > analysis_time_max) continue;

                double swap_dt = swaprecords[nf][nsj].time - swaprecords[nf][nsi].time;
                double swap_dr = (swaprecords[nf][nsj].bondCenter - swaprecords[nf][nsi].bondCenter).length();
                if (swap_dt > 0.00)
                {
                    swapStat_ijCount_whistdt.add(log10(swap_dt), 1.00);
                    int rbin = swap_dr / swapStat_ijCount_rbinsize;
                    if (rbin >= 0 && rbin < swapStat_ij_dr_nbins)
                        swapStat_ijCount_whistdtdr[rbin].add(log10(swap_dt), 1.00);

                    double swap_dtt = swap_dt / swaprecords[nf][nsj].time;
                    hist_ij_logdtdr.add(log10(swap_dt), swap_dr);
                    if (swap_dr > 0.00) hist_ij_logdtlogdr.add(log10(swap_dt), log10(swap_dr));
                }
                double dtheta = log10(swaprecords[nf][nsj].time) - log10(swaprecords[nf][nsi].time);
                if (dtheta > 0.00)
                {
                    swapStat_ijCount_whistdtheta.add(log10(dtheta), 1.00);
                    int rbin = swap_dr / swapStat_ijCount_rbinsize;
                    if (rbin >= 0 && rbin < swapStat_ij_dr_nbins)
                        swapStat_ijCount_whistdthetadr[rbin].add(log10(dtheta), 1.00);
                }
                swapStat_ijCount_whistdr.add(swap_dr, 1.00);
            }
            for (int nsj = 0; nsj < swaprecords[nf].size(); nsj++)
            {
                if (nsi == nsj) continue;
                if (swaprecords[nf][nsj].time < analysis_time_min) continue;
                if (analysis_time_max > 0.00 && swaprecords[nf][nsj].time > analysis_time_max) continue;

                double swap_dr = (swaprecords[nf][nsj].bondCenter - swaprecords[nf][nsi].bondCenter).length();
                swapStat_ijCount_whistdr_alltime.add(swap_dr, 1.00);
            }
        }
    }

    // INTRA SUBAVALANCHE
    cout << "Intra-subavalanche swaps\n";
    vector<v3_math::Histogramm<double>> intrasubav_swap_logdt_epoch(swapStat_nepoch * swapStat_nepochfrac + 1);
    vector<v3_math::Histogramm<double>> intrasubav_swap_dtheta_epoch(swapStat_nepoch * swapStat_nepochfrac + 1);
    vector<v3_math::CHistogramm<double>> intrasubav_swap_cdtheta_epoch(swapStat_nepoch * swapStat_nepochfrac + 1);
    vector<v3_math::CHistogramm<double>> intrasubav_swap_cdtt_epoch(swapStat_nepoch * swapStat_nepochfrac + 1);
    v3_math::Histogramm<double> intrasubav_swap_histdr(swapStat_dr_min, swapStat_dr_max, swapStat_dr_nbins, v3_math::HISTMODE_CUTOFF);
    v3_math::Histogramm<double> intrasubav_swap_histlogdr(swapStat_logdr_min, swapStat_logdr_max, swapStat_logdr_nbins, v3_math::HISTMODE_CUTOFF);
    v3_math::Histogramm<double> intrasubav_swap_histdrclosest(swapStat_dr_min, swapStat_dr_max, swapStat_dr_nbins, v3_math::HISTMODE_CUTOFF);
    v3_math::Histogramm<double> intrasubav_swap_histlogdrclosest(swapStat_logdr_min, swapStat_logdr_max, swapStat_logdr_nbins, v3_math::HISTMODE_CUTOFF);
    v3_math::WHistogramm<double> intrasubav_swap_ijCount_whistdr(swapStat_ij_dr_min, swapStat_ij_dr_max, swapStat_ij_dr_nbins, v3_math::HISTMODE_EXT),
        intrasubav_swap_ijNoweight_whistdr(swapStat_ij_dr_min, swapStat_ij_dr_max, swapStat_ij_dr_nbins, v3_math::HISTMODE_EXT);
    v3_math::WHistogramm<double> intrasubav_swap_ijCount_whistdr_alltime(swapStat_ij_dr_min, swapStat_ij_dr_max, swapStat_ij_dr_nbins, v3_math::HISTMODE_EXT),
        intrasubav_swap_ijNoweight_whistdr_alltime(swapStat_ij_dr_min, swapStat_ij_dr_max, swapStat_ij_dr_nbins, v3_math::HISTMODE_EXT);
    for (int q = 0; q < swapStat_nepoch * swapStat_nepochfrac + 1; q++)
    {
        intrasubav_swap_logdt_epoch[q].change(-3.00, 7.00, 10 * 3, v3_math::HISTMODE_EXT);
        intrasubav_swap_dtheta_epoch[q].change(0.00, 0.50, 100, v3_math::HISTMODE_EXT);
        intrasubav_swap_cdtheta_epoch[q].change(dthetaBinBorders, v3_math::CHISTMODE_EXT);
        intrasubav_swap_cdtt_epoch[q].change(dthetaBinBorders, v3_math::CHISTMODE_EXT);
    }
    vector<double> fractDimCubeLengths;
    fractDimCubeLengths.push_back(1000.0);
    fractDimCubeLengths.push_back(900.0);
    fractDimCubeLengths.push_back(800.0);
    fractDimCubeLengths.push_back(700.0);
    fractDimCubeLengths.push_back(600.0);
    fractDimCubeLengths.push_back(500.0);
    fractDimCubeLengths.push_back(450.0);
    fractDimCubeLengths.push_back(400.0);
    fractDimCubeLengths.push_back(350.0);
    fractDimCubeLengths.push_back(300.0);
    fractDimCubeLengths.push_back(250.0);
    fractDimCubeLengths.push_back(200.0);
    fractDimCubeLengths.push_back(175.0);
    fractDimCubeLengths.push_back(150.0);
    fractDimCubeLengths.push_back(125.0);
    fractDimCubeLengths.push_back(100.0);
    fractDimCubeLengths.push_back(90.0);
    fractDimCubeLengths.push_back(80.0);
    fractDimCubeLengths.push_back(70.0);
    fractDimCubeLengths.push_back(60.0);
    fractDimCubeLengths.push_back(50.0);
    fractDimCubeLengths.push_back(40.0);
    fractDimCubeLengths.push_back(30.0);
    fractDimCubeLengths.push_back(25.0);
    fractDimCubeLengths.push_back(20.0);
    fractDimCubeLengths.push_back(15.0);
    fractDimCubeLengths.push_back(12.5);
    fractDimCubeLengths.push_back(10.0);
    fractDimCubeLengths.push_back(7.5);
    fractDimCubeLengths.push_back(5.0);
    int nlengths = fractDimCubeLengths.size() < subav_fractDim_nLengths ? fractDimCubeLengths.size() : subav_fractDim_nLengths;
    vector<vector<vector<double>>> subav_FractDim_N(nfiles), subav_FractDim_S(nfiles), subav_FractDim_L(nfiles);
    vector<v3_math::AvgErr<double>> subav_FractDim_N_avg(nlengths);
    vector<v3_math::AvgErr<double>> subav_FractDim_NNmax_avg(nlengths);
    vector<v3_math::AvgErr<double>> subav_FractDim_logN_avg(nlengths);
    vector<v3_math::AvgErr<double>> subav_FractDim_logS_avg(nlengths);
    vector<v3_math::AvgErr<double>> subav_FractDim_logL_avg(nlengths);
    vector<v3_math::AvgErr<double>> subav_FractDim_logStoL_avg(nlengths);
    vector<v3_math::Histogramm<double>> subav_FractDim_logS_hist(nlengths);
    vector<v3_math::Histogramm<double>> subav_FractDim_logL_hist(nlengths);
    vector<v3_math::Histogramm<double>> subav_FractDim_logStoL_hist(nlengths);
    for (int q = 0; q < nlengths; q++)
    {
        subav_FractDim_logS_hist[q].change(subav_fractDim_logS_min, subav_fractDim_logS_max, subav_fractDim_logS_nbins);
        subav_FractDim_logL_hist[q].change(subav_fractDim_logL_min, subav_fractDim_logL_max, subav_fractDim_logL_nbins);
        subav_FractDim_logStoL_hist[q].change(subav_fractDim_logStoL_min, subav_fractDim_logStoL_max, subav_fractDim_logStoL_nbins);
    }
    for (int nf = 0; nf < nfiles; nf++)
    {
        cout << "\nFile " << nf + 1 << " / " << nfiles<<"\n";
        cout << "boxLimits: x in " << boxLimits[nf].L.x << " to " << boxLimits[nf].R.x << ", y in " << boxLimits[nf].L.y << " to " << boxLimits[nf].R.y << "\n";

        for(int nsa = 0; nsa < subavalanches[nf].size(); nsa++)
        {
            //cout << "\nsubav "<<it->first;
            //cout.flush();
            for (int nss = 1; nss < subavalanches[nf][nsa].swaprecords_time_ID.size(); nss++)
            {
                int swapID = subavalanches[nf][nsa].swaprecords_time_ID[nss].second;
                int prevswapID = subavalanches[nf][nsa].swaprecords_time_ID[nss - 1].second;

                if (swaprecords[nf][swapID].time < analysis_time_min) continue;
                if (analysis_time_max > 0.00 && swaprecords[nf][swapID].time > analysis_time_max) continue;

                double swapTime_log10 = log10(swaprecords[nf][swapID].time);
                double dtime = swaprecords[nf][swapID].time - swaprecords[nf][prevswapID].time;
                double dtheta = swapTime_log10 - log10(swaprecords[nf][prevswapID].time);
                double swap_dr = (swaprecords[nf][swapID].bondCenter - swaprecords[nf][prevswapID].bondCenter).length();
                intrasubav_swap_histdr.add(swap_dr);
                if (swap_dr > 0.00) intrasubav_swap_histlogdr.add(log10(swap_dr));

                int ne = swapTime_log10 - swapStat_epochStart_log10;
                int nef = (swapTime_log10 - swapStat_epochStart_log10) * double(swapStat_nepochfrac);

                if (dtime > 0.00)
                {
                    double log10dt = log10(dtime);
                    intrasubav_swap_logdt_epoch[0].add(log10dt);
                    intrasubav_swap_dtheta_epoch[0].add(dtheta);
                    intrasubav_swap_cdtheta_epoch[0].add(dtheta);
                    if (nef >= 0 && nef < swapStat_nepoch * swapStat_nepochfrac)
                    {
                        intrasubav_swap_logdt_epoch[nef + 1].add(log10dt);
                        intrasubav_swap_dtheta_epoch[nef + 1].add(dtheta);
                        intrasubav_swap_cdtheta_epoch[nef + 1].add(dtheta);
                    }
                }
                if (swaprecords[nf][swapID].time > 0.00)
                {
                    double dtt = dtime / swaprecords[nf][swapID].time;
                    intrasubav_swap_cdtt_epoch[0].add(dtt);
                    if (nef >= 0 && nef < swapStat_nepoch * swapStat_nepochfrac)
                    {
                        intrasubav_swap_cdtt_epoch[nef + 1].add(dtt);
                    }
                }
            }
            for (int nss = 0; nss < subavalanches[nf][nsa].swaprecords_time_ID.size(); nss++)
            {
                int swapID = subavalanches[nf][nsa].swaprecords_time_ID[nss].second;

                if (swaprecords[nf][swapID].time < analysis_time_min) continue;
                if (analysis_time_max > 0.00 && swaprecords[nf][swapID].time > analysis_time_max) continue;

                v3_math::planar::point<double> boxLeft(system_posLeft_x, system_posLeft_y), boxRight(system_posLeft_x + swaprecords[nf][swapID].xsize, system_posLeft_y + swaprecords[nf][swapID].ysize);
                v3_math::planar::rectangle<double> box(boxLeft, boxRight);
                v3_math::planar::point<double> bondCenter(swaprecords[nf][swapID].bondCenter.x, swaprecords[nf][swapID].bondCenter.y);
                for (int rbin = 0; rbin < swapStat_ij_dr_nbins; rbin++)
                {
                    double rin = swapStat_ij_dr_min + double(rbin) * swapStat_ijCount_rbinsize;
                    double rout = swapStat_ij_dr_min + double(rbin + 1) * swapStat_ijCount_rbinsize;
                    v3_math::planar::circle<double> Cin(bondCenter, rin), Cout(bondCenter, rout);
                    double Sin = v3_math::planar::intersectarea(box, Cin);
                    double Sout = v3_math::planar::intersectarea(box, Cout);
                    double Sbin = Sout - Sin;
                    intrasubav_swap_ijNoweight_whistdr.add((rin + rout) * 0.50, Sbin);
                    intrasubav_swap_ijNoweight_whistdr_alltime.add((rin + rout) * 0.50, Sbin);
                }
                
                for (int nssj = nss + 1; nssj < subavalanches[nf][nsa].swaprecords_time_ID.size(); nssj++)
                {
                    int nextswapID = subavalanches[nf][nsa].swaprecords_time_ID[nssj].second;
                    double swap_dr = (swaprecords[nf][swapID].bondCenter - swaprecords[nf][nextswapID].bondCenter).length();
                    intrasubav_swap_ijCount_whistdr.add(swap_dr, 1.00);
                }
                for (int nssj = 0; nssj < subavalanches[nf][nsa].swaprecords_time_ID.size(); nssj++)
                {
                    if (nssj == nss) continue;
                    int nextswapID = subavalanches[nf][nsa].swaprecords_time_ID[nssj].second;
                    double swap_dr = (swaprecords[nf][swapID].bondCenter - swaprecords[nf][nextswapID].bondCenter).length();
                    intrasubav_swap_ijCount_whistdr_alltime.add(swap_dr, 1.00);
                }
            }
            for (int nss = 0; nss < subavalanches[nf][nsa].swaprecords_time_ID.size(); nss++)
            {
                int swapID = subavalanches[nf][nsa].swaprecords_time_ID[nss].second;
                double dr_closest = -1.00;
                for (int q = 0; q < subavalanches[nf][nsa].swaprecords_time_ID.size(); q++)
                {
                    if (nss == q) continue;
                    int qswapID = subavalanches[nf][nsa].swaprecords_time_ID[q].second;
                    double dt = swaprecords[nf][swapID].time - swaprecords[nf][qswapID].time;
                    if (dt < 0.00) dt *= -1.00;
                    if (dt > subavalanche_dt_threshold) continue;
                    double swap_dr = (swaprecords[nf][swapID].bondCenter - swaprecords[nf][qswapID].bondCenter).length();
                    if (dr_closest<0.00 || dr_closest>swap_dr) dr_closest = swap_dr;
                }
                intrasubav_swap_histdrclosest.add(dr_closest);
                if (dr_closest > 0.00) intrasubav_swap_histlogdrclosest.add(log10(dr_closest));
            }         
        } // for (map<int, vector<int>>::const_iterator it = subavalanche_swaps[nf].begin(); it != subavalanche_swaps[nf].end(); it++)
        
        // Fractal dimension
        cout << "FD\n";
        subav_FractDim_N[nf].resize(subavalanches[nf].size());
        subav_FractDim_S[nf].resize(subavalanches[nf].size());
        subav_FractDim_L[nf].resize(subavalanches[nf].size());
        for (int nsa = 0; nsa < subavalanches[nf].size(); nsa++)
        {
            subav_FractDim_N[nf][nsa].resize(nlengths);
            subav_FractDim_S[nf][nsa].resize(nlengths);
            subav_FractDim_L[nf][nsa].resize(nlengths);
        }
        for (int fdn = 0; fdn < nlengths; fdn++)
        {
            double fdl = fractDimCubeLengths[fdn];
            double nx = (boxLimits[nf].R.x - boxLimits[nf].L.x) / fdl;
            double ny = (boxLimits[nf].R.y - boxLimits[nf].L.y) / fdl;
            int cnx = nx, cny = ny;
            if (nx > double(cnx)) cnx++;
            if (ny > double(cny)) cny++;
            cnx++;
            cny++;
            v3_math::Matrix<bool> cubes(cnx, cny);
            cout << "FDN = " << fdn << ", FDL = " << fdl << ", cubes " << cnx << " x " << cny << "\n";

            for (int nsa = 0; nsa < subavalanches[nf].size(); nsa++)
            {
                //if (nsa != 1294) continue;

                // Check general time limit
                if (subavalanches[nf][nsa].subavTimeP50 < analysis_time_min) continue;
                if (analysis_time_max > 0.00 && subavalanches[nf][nsa].subavTimeP50 > analysis_time_max) continue;

                //cout << nsa << "\n";
                cubes.set(false);
                for (int nss = 0; nss < subavalanches[nf][nsa].swaprecords_time_ID.size(); nss++)
                {
                    int swapID = subavalanches[nf][nsa].swaprecords_time_ID[nss].second;
                    int snx = (swaprecords[nf][swapID].bondCenter.x - boxLimits[nf].L.x) / fdl;
                    int sny = (swaprecords[nf][swapID].bondCenter.y - boxLimits[nf].L.y) / fdl;
                    cubes(snx, sny) = true;
                }
                int ncubes = 0, nsides = 0;
                for (int ix = 0; ix < cnx; ix++)
                {
                    for (int iy = 0; iy < cny; iy++)
                    {
                        if (cubes(ix, iy) == true) ncubes++;
                        if (ix == 0) nsides++;
                        else if (cubes(ix - 1, iy) == false) nsides++;
                        if (ix == cnx - 1) nsides++;
                        else if (cubes(ix + 1, iy) == false) nsides++;
                        if (iy == 0) nsides++;
                        else if (cubes(ix, iy - 1) == false) nsides++;
                        if (iy == cny - 1) nsides++;
                        else if (cubes(ix, iy + 1) == false) nsides++;
                    }
                }
                
                subav_FractDim_N[nf][nsa][fdn] = ncubes;
                subav_FractDim_S[nf][nsa][fdn] = double(ncubes) * fdl * fdl;
                subav_FractDim_L[nf][nsa][fdn] = double(nsides) * fdl;
                if (subavalanches[nf][nsa].nswaps >= subav_fractDim_nswapsmin && subavalanches[nf][nsa].nswaps <= subav_fractDim_nswapsmax &&
                    subavalanches[nf][nsa].subavRg >= subav_fractDim_Rgmin && subavalanches[nf][nsa].subavRg <= subav_fractDim_Rgmax)
                {
                    subav_FractDim_N_avg[fdn].add(ncubes);
                    subav_FractDim_logN_avg[fdn].add(log10(ncubes));
                    subav_FractDim_logS_avg[fdn].add(log10(subav_FractDim_S[nf][nsa][fdn]));
                    subav_FractDim_logL_avg[fdn].add(log10(subav_FractDim_L[nf][nsa][fdn]));
                    subav_FractDim_logStoL_avg[fdn].add(log10(subav_FractDim_S[nf][nsa][fdn] / subav_FractDim_L[nf][nsa][fdn]));

                    subav_FractDim_logS_hist[fdn].add(log10(subav_FractDim_S[nf][nsa][fdn]));
                    subav_FractDim_logL_hist[fdn].add(log10(subav_FractDim_L[nf][nsa][fdn]));
                    subav_FractDim_logStoL_hist[fdn].add(log10(subav_FractDim_S[nf][nsa][fdn] / subav_FractDim_L[nf][nsa][fdn]));
                }
            }
        }
        for (int fdn = 0; fdn < nlengths; fdn++)
        {
            for (int nsa = 0; nsa < subavalanches[nf].size(); nsa++)
            {
                if (subavalanches[nf][nsa].subavTimeP50 < analysis_time_min) continue;
                if (analysis_time_max > 0.00 && subavalanches[nf][nsa].subavTimeP50 > analysis_time_max) continue;
                if (subavalanches[nf][nsa].nswaps < subav_fractDim_nswapsmin || subavalanches[nf][nsa].nswaps > subav_fractDim_nswapsmax ||
                    subavalanches[nf][nsa].subavRg < subav_fractDim_Rgmin || subavalanches[nf][nsa].subavRg > subav_fractDim_Rgmax) continue;
                
                subav_FractDim_NNmax_avg[fdn].add(double(subav_FractDim_N[nf][nsa][fdn]) / double(subav_FractDim_N[nf][nsa][nlengths - 1]));
            }
        }
        cout << "subavs+ ";
        cout.flush();
    }

    // INTRA AVALANCHE
    cout << "Intra-avalanche swaps\n";
    v3_math::Histogramm<double> intraav_swap_histdr(swapStat_dr_min, swapStat_dr_max, swapStat_dr_nbins, v3_math::HISTMODE_CUTOFF);
    v3_math::Histogramm<double> intraav_swap_histlogdr(log10(swapStat_logdr_min), log10(swapStat_logdr_max), swapStat_logdr_nbins, v3_math::HISTMODE_CUTOFF);
    v3_math::WHistogramm<double> intraav_swap_ijCount_whistdr(swapStat_ij_dr_min, swapStat_ij_dr_max, swapStat_ij_dr_nbins, v3_math::HISTMODE_EXT),
        intraav_swap_ijNoweight_whistdr(swapStat_ij_dr_min, swapStat_ij_dr_max, swapStat_ij_dr_nbins, v3_math::HISTMODE_EXT);
    v3_math::WHistogramm<double> intraav_swap_ijCount_whistdr_alltime(swapStat_ij_dr_min, swapStat_ij_dr_max, swapStat_ij_dr_nbins, v3_math::HISTMODE_EXT),
        intraav_swap_ijNoweight_whistdr_alltime(swapStat_ij_dr_min, swapStat_ij_dr_max, swapStat_ij_dr_nbins, v3_math::HISTMODE_EXT);
    vector<vector<vector<double>>> av_FractDim_N(nfiles), av_FractDim_S(nfiles), av_FractDim_L(nfiles);
    vector<v3_math::AvgErr<double>> av_FractDim_N_avg(nlengths);
    vector<v3_math::AvgErr<double>> av_FractDim_logN_avg(nlengths);
    vector<v3_math::AvgErr<double>> av_FractDim_logS_avg(nlengths);
    vector<v3_math::AvgErr<double>> av_FractDim_logL_avg(nlengths);
    vector<v3_math::AvgErr<double>> av_FractDim_logStoL_avg(nlengths);
    vector<v3_math::Histogramm<double>> av_FractDim_logS_hist(nlengths);
    vector<v3_math::Histogramm<double>> av_FractDim_logL_hist(nlengths);
    vector<v3_math::Histogramm<double>> av_FractDim_logStoL_hist(nlengths);
    vector<vector<vector<double>>> av_subavFractDim_N(nfiles), av_subavFractDim_S(nfiles), av_subavFractDim_L(nfiles);
    vector<v3_math::AvgErr<double>> av_subavFractDim_N_avg(nlengths);
    vector<v3_math::AvgErr<double>> av_subavFractDim_logN_avg(nlengths);
    vector<v3_math::AvgErr<double>> av_subavFractDim_logS_avg(nlengths);
    vector<v3_math::AvgErr<double>> av_subavFractDim_logL_avg(nlengths);
    vector<v3_math::AvgErr<double>> av_subavFractDim_logStoL_avg(nlengths);
    vector<v3_math::Histogramm<double>> av_subavFractDim_logS_hist(nlengths);
    vector<v3_math::Histogramm<double>> av_subavFractDim_logL_hist(nlengths);
    vector<v3_math::Histogramm<double>> av_subavFractDim_logStoL_hist(nlengths);
    for (int q = 0; q < nlengths; q++)
    {
        av_FractDim_logS_hist[q].change(subav_fractDim_logS_min, subav_fractDim_logS_max, subav_fractDim_logS_nbins);
        av_FractDim_logL_hist[q].change(subav_fractDim_logL_min, subav_fractDim_logL_max, subav_fractDim_logL_nbins);
        av_FractDim_logStoL_hist[q].change(subav_fractDim_logStoL_min, subav_fractDim_logStoL_max, subav_fractDim_logStoL_nbins);
        av_subavFractDim_logS_hist[q].change(subav_fractDim_logS_min, subav_fractDim_logS_max, subav_fractDim_logS_nbins);
        av_subavFractDim_logL_hist[q].change(subav_fractDim_logL_min, subav_fractDim_logL_max, subav_fractDim_logL_nbins);
        av_subavFractDim_logStoL_hist[q].change(subav_fractDim_logStoL_min, subav_fractDim_logStoL_max, subav_fractDim_logStoL_nbins);
    }
    for (int nf = 0; nf < nfiles; nf++)
    {
        cout << "\nFile " << nf + 1 << " / " << nfiles;
        cout.flush();

        for (int nav = 0; nav < avalanches[nf].size(); nav++)
        {
            for (int ns = 1; ns < avalanches[nf][nav].swaprecords_time_ID.size(); ns++)
            {
                int swapID = avalanches[nf][nav].swaprecords_time_ID[ns].second;

                if (swaprecords[nf][swapID].time < analysis_time_min) continue;
                if (analysis_time_max > 0.00 && swaprecords[nf][swapID].time > analysis_time_max) continue;

                int prevswapID = avalanches[nf][nav].swaprecords_time_ID[ns - 1].second;

                double swap_dr = (swaprecords[nf][swapID].bondCenter - swaprecords[nf][prevswapID].bondCenter).length();
                intraav_swap_histdr.add(swap_dr);
                if (swap_dr > 0.00) intraav_swap_histlogdr.add(log10(swap_dr));
            }
            for (int nsi = 0; nsi < avalanches[nf][nav].swaprecords_time_ID.size(); nsi++)
            {
                int swapID = avalanches[nf][nav].swaprecords_time_ID[nsi].second;

                if (swaprecords[nf][swapID].time < analysis_time_min) continue;
                if (analysis_time_max > 0.00 && swaprecords[nf][swapID].time > analysis_time_max) continue;

                v3_math::planar::point<double> boxLeft(system_posLeft_x, system_posLeft_y), boxRight(system_posLeft_x + swaprecords[nf][swapID].xsize, system_posLeft_y + swaprecords[nf][swapID].ysize);
                v3_math::planar::rectangle<double> box(boxLeft, boxRight);
                v3_math::planar::point<double> bondCenter(swaprecords[nf][swapID].bondCenter.x, swaprecords[nf][swapID].bondCenter.y);
                for (int rbin = 0; rbin < swapStat_ij_dr_nbins; rbin++)
                {
                    double rin = swapStat_ij_dr_min + double(rbin) * swapStat_ijCount_rbinsize;
                    double rout = swapStat_ij_dr_min + double(rbin + 1) * swapStat_ijCount_rbinsize;
                    v3_math::planar::circle<double> Cin(bondCenter, rin), Cout(bondCenter, rout);
                    double Sin = v3_math::planar::intersectarea(box, Cin);
                    double Sout = v3_math::planar::intersectarea(box, Cout);
                    double Sbin = Sout - Sin;
                    intraav_swap_ijNoweight_whistdr.add((rin + rout) * 0.50, Sbin);
                    intraav_swap_ijNoweight_whistdr_alltime.add((rin + rout) * 0.50, Sbin);
                }
                for (int nsj = nsi + 1; nsj < avalanches[nf][nav].swaprecords_time_ID.size(); nsj++)
                {
                    int nextswapID = avalanches[nf][nav].swaprecords_time_ID[nsj].second;
                    double swap_dr = (swaprecords[nf][swapID].bondCenter - swaprecords[nf][nextswapID].bondCenter).length();
                    intraav_swap_ijCount_whistdr.add(swap_dr, 1.00);
                }
                for (int nsj = 0; nsj < avalanches[nf][nav].swaprecords_time_ID.size(); nsj++)
                {
                    if (nsj == nsi) continue;
                    int nextswapID = avalanches[nf][nav].swaprecords_time_ID[nsj].second;
                    double swap_dr = (swaprecords[nf][swapID].bondCenter - swaprecords[nf][nextswapID].bondCenter).length();
                    intraav_swap_ijCount_whistdr_alltime.add(swap_dr, 1.00);
                }
            }
        }

        // Fractal dimension
        cout << "FD\n";
        av_FractDim_N[nf].resize(avalanches[nf].size());
        av_FractDim_S[nf].resize(avalanches[nf].size());
        av_FractDim_L[nf].resize(avalanches[nf].size());
        av_subavFractDim_N[nf].resize(avalanches[nf].size());
        av_subavFractDim_S[nf].resize(avalanches[nf].size());
        av_subavFractDim_L[nf].resize(avalanches[nf].size());
        for (int nav = 0; nav < avalanches[nf].size(); nav++)
        {
            av_FractDim_N[nf][nav].resize(nlengths);
            av_FractDim_S[nf][nav].resize(nlengths);
            av_FractDim_L[nf][nav].resize(nlengths);
            av_subavFractDim_N[nf][nav].resize(nlengths);
            av_subavFractDim_S[nf][nav].resize(nlengths);
            av_subavFractDim_L[nf][nav].resize(nlengths);
        }
        for (int fdn = 0; fdn < nlengths; fdn++)
        {
            double fdl = fractDimCubeLengths[fdn];
            double nx = (boxLimits[nf].R.x - boxLimits[nf].L.x) / fdl;
            double ny = (boxLimits[nf].R.y - boxLimits[nf].L.y) / fdl;
            int cnx = nx, cny = ny;
            if (nx > double(cnx)) cnx++;
            if (ny > double(cny)) cny++;
            cnx++;
            cny++;
            v3_math::Matrix<bool> cubes(cnx, cny);
            cout << "FDN = " << fdn << ", FDL = " << fdl << ", cubes " << cnx << " x " << cny << "\n";

            for (int nav = 0; nav < avalanches[nf].size(); nav++)
            {
                //if (nav != 732) continue;

                if (avalanches[nf][nav].avTimeP50 < analysis_time_min) continue;
                if (analysis_time_max > 0.00 && avalanches[nf][nav].avTimeP50 > analysis_time_max) continue;

                //cout << nsa << "\n";
                cubes.set(false);
                for (int nsi = 0; nsi < avalanches[nf][nav].swaprecords_time_ID.size(); nsi++)
                {
                    int swapID = avalanches[nf][nav].swaprecords_time_ID[nsi].second;

                    int snx = (swaprecords[nf][swapID].bondCenter.x - boxLimits[nf].L.x) / fdl;
                    int sny = (swaprecords[nf][swapID].bondCenter.y - boxLimits[nf].L.y) / fdl;
                    cubes(snx, sny) = true;
                }
                int ncubes = 0, nsides = 0;
                for (int ix = 0; ix < cnx; ix++)
                {
                    for (int iy = 0; iy < cny; iy++)
                    {
                        if (cubes(ix, iy) == true) ncubes++;
                        if (ix == 0) nsides++;
                        else if (cubes(ix - 1, iy) == false) nsides++;
                        if (ix == cnx - 1) nsides++;
                        else if (cubes(ix + 1, iy) == false) nsides++;
                        if (iy == 0) nsides++;
                        else if (cubes(ix, iy - 1) == false) nsides++;
                        if (iy == cny - 1) nsides++;
                        else if (cubes(ix, iy + 1) == false) nsides++;
                    }
                }

                av_FractDim_N[nf][nav][fdn] = ncubes;
                av_FractDim_S[nf][nav][fdn] = double(ncubes) * fdl * fdl;
                av_FractDim_L[nf][nav][fdn] = double(nsides) * fdl;
                if (avalanches[nf][nav].avNswaps >= av_fractDim_nswapsmin && avalanches[nf][nav].avNswaps <= av_fractDim_nswapsmax)
                {
                    av_FractDim_N_avg[fdn].add(ncubes);
                    av_FractDim_logN_avg[fdn].add(log10(ncubes));
                    av_FractDim_logS_avg[fdn].add(log10(av_FractDim_S[nf][nav][fdn]));
                    av_FractDim_logL_avg[fdn].add(log10(av_FractDim_L[nf][nav][fdn]));
                    av_FractDim_logStoL_avg[fdn].add(log10(av_FractDim_S[nf][nav][fdn] / av_FractDim_L[nf][nav][fdn]));

                    av_FractDim_logS_hist[fdn].add(log10(av_FractDim_S[nf][nav][fdn]));
                    av_FractDim_logL_hist[fdn].add(log10(av_FractDim_L[nf][nav][fdn]));
                    av_FractDim_logStoL_hist[fdn].add(log10(av_FractDim_S[nf][nav][fdn] / av_FractDim_L[nf][nav][fdn]));
                }

                cubes.set(false);
                for (int nsa = 0; nsa < avalanches[nf][nav].subavalanches_timeP50_ID.size(); nsa++)
                {
                    int subavID = avalanches[nf][nav].subavalanches_timeP50_ID[nsa].second;

                    int snx = (subavalanches[nf][subavID].subavCenterX - boxLimits[nf].L.x) / fdl;
                    int sny = (subavalanches[nf][subavID].subavCenterY - boxLimits[nf].L.y) / fdl;
                    cubes(snx, sny) = true;
                }
                ncubes = 0;
                nsides = 0;
                for (int ix = 0; ix < cnx; ix++)
                {
                    for (int iy = 0; iy < cny; iy++)
                    {
                        if (cubes(ix, iy) == true) ncubes++;
                        if (ix == 0) nsides++;
                        else if (cubes(ix - 1, iy) == false) nsides++;
                        if (ix == cnx - 1) nsides++;
                        else if (cubes(ix + 1, iy) == false) nsides++;
                        if (iy == 0) nsides++;
                        else if (cubes(ix, iy - 1) == false) nsides++;
                        if (iy == cny - 1) nsides++;
                        else if (cubes(ix, iy + 1) == false) nsides++;
                    }
                }

                av_subavFractDim_N[nf][nav][fdn] = ncubes;
                av_subavFractDim_S[nf][nav][fdn] = double(ncubes) * fdl * fdl;
                av_subavFractDim_L[nf][nav][fdn] = double(nsides) * fdl;
                if (avalanches[nf][nav].avNswaps >= av_fractDim_nswapsmin && avalanches[nf][nav].avNswaps <= av_fractDim_nswapsmax)
                {
                    av_subavFractDim_N_avg[fdn].add(ncubes);
                    av_subavFractDim_logN_avg[fdn].add(log10(ncubes));
                    av_subavFractDim_logS_avg[fdn].add(log10(av_FractDim_S[nf][nav][fdn]));
                    av_subavFractDim_logL_avg[fdn].add(log10(av_FractDim_L[nf][nav][fdn]));
                    av_subavFractDim_logStoL_avg[fdn].add(log10(av_subavFractDim_S[nf][nav][fdn] / av_subavFractDim_L[nf][nav][fdn]));

                    av_subavFractDim_logS_hist[fdn].add(log10(av_subavFractDim_S[nf][nav][fdn]));
                    av_subavFractDim_logL_hist[fdn].add(log10(av_subavFractDim_L[nf][nav][fdn]));
                    av_subavFractDim_logStoL_hist[fdn].add(log10(av_subavFractDim_S[nf][nav][fdn] / av_subavFractDim_L[nf][nav][fdn]));
                }
            }
        }
        cout << "avs+ ";
        cout.flush();
    }


    cout << "done.\n";
    /** / Collect swap statistics **/

/**
* This block calculates statistics on subavalanches/cascades
**/
    /** $ Collect subavalanche statistics **/
    cout << "\nCollect subavalanche statistics...";
    // ALL SUBAVALANCHES
    cout << "All subavalanches\n";
    // Single subavalanche statistics
    v3_math::Histogramm<double> subavStat_histNswaps(subavStat_nswaps_log10_min, subavStat_nswaps_log10_max, subavStat_nswaps_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm<double> subavStat_histDuration(subavStat_duration_log10_min, subavStat_duration_log10_max, subavStat_duration_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::WHistogramm<double> subavStat_whistDuration(subavStat_duration_log10_min, subavStat_duration_log10_max, subavStat_duration_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm<double> subavStat_histSwapRate(subavStat_swaprate_log10_min, subavStat_swaprate_log10_max, subavStat_swaprate_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm<double> subavStat_histDiameter(subavStat_diameter_log10_min, subavStat_diameter_log10_max, subavStat_diameter_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm<double> subavStat_histRg(subavStat_rg_log10_min, subavStat_rg_log10_max, subavStat_rg_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm<double> subavStat_histI12(subavStat_I12_log10_min, subavStat_I12_log10_max, subavStat_I12_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::WHistogramm<double> subavStat_whistI12(subavStat_I12_log10_min, subavStat_I12_log10_max, subavStat_I12_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm2D<double> subavStat_histNswapsDiameter(subavStat_nswaps_log10_min, subavStat_nswaps_log10_max, subavStat_diameter_log10_min, 
        subavStat_diameter_log10_max, subavStat_nswaps_log10_nbins, subavStat_diameter_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm2D<double> subavStat_histNswapsRg(subavStat_nswaps_log10_min, subavStat_nswaps_log10_max, subavStat_rg_log10_min, subavStat_rg_log10_max,
        subavStat_nswaps_log10_nbins, subavStat_rg_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm2D<double> subavStat_histNswapsDuration(subavStat_nswaps_log10_min, subavStat_nswaps_log10_max, 
        subavStat_duration_log10_min, subavStat_duration_log10_max,
        subavStat_nswaps_log10_nbins, subavStat_duration_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm2D<double> subavStat_histNswapsDirDiameter(subavStat_nswaps_log10_min, subavStat_nswaps_log10_max, 
        subavStat_diameter_log10_min, subavStat_diameter_log10_max, subavStat_nswaps_log10_nbins, subavStat_diameter_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm2D<double> subavStat_histNswapsDirRg(subavStat_nswaps_log10_min, subavStat_nswaps_log10_max, subavStat_rg_log10_min, subavStat_rg_log10_max,
        subavStat_nswaps_log10_nbins, subavStat_rg_log10_nbins, v3_math::HISTMODE_EXT);
    v3_phys::path_analyser subav_path_analyzer(subav_path_setup);
    // Inter-subavalanche statistics
    v3_math::Histogramm<double> subavStat_histdtP50(subavStat_dtP50_log10_min, subavStat_dtP50_log10_max, subavStat_dtP50_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm<double> subavStat_histdttP50(subavStat_dttP50_log10_min, subavStat_dttP50_log10_max, subavStat_dttP50_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::WHistogramm<double> subavStat_histdttP50_weighted(subavStat_dttP50_log10_min, subavStat_dttP50_log10_max, subavStat_dttP50_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm<double> subavStat_histdtETS(subavStat_dtP50_log10_min, subavStat_dtP50_log10_max, subavStat_dtP50_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm<double> subavStat_histdttETS(subavStat_dttP50_log10_min, subavStat_dttP50_log10_max, subavStat_dttP50_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm<double> subavStat_histdr(subavStat_dr_min, subavStat_dr_max, subavStat_dr_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm<double> subavStat_histlogdr(subavStat_dr_log10_min, subavStat_dr_log10_max, subavStat_dr_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm2D<double> subavStat_histdtP50dr(subavStat_dtP50_log10_min, subavStat_dtP50_log10_max, subavStat_dr_min, subavStat_dr_max,
        subavStat_dtP50_log10_nbins, subavStat_dr_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm2D<double> subavStat_histdttP50dr(subavStat_dttP50_log10_min, subavStat_dttP50_log10_max, subavStat_dr_min, subavStat_dr_max,
        subavStat_dttP50_log10_nbins, subavStat_dr_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm2D<double> subavStat_histdtP50logdr(subavStat_dtP50_log10_min, subavStat_dtP50_log10_max, subavStat_dr_log10_min, subavStat_dr_log10_max,
        subavStat_dtP50_log10_nbins, subavStat_dr_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm2D<double> subavStat_histdttP50logdr(subavStat_dttP50_log10_min, subavStat_dttP50_log10_max, subavStat_dr_log10_min, subavStat_dr_log10_max,
        subavStat_dttP50_log10_nbins, subavStat_dr_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm2D<double> subavStat_histdtETSdr(subavStat_dtP50_log10_min, subavStat_dtP50_log10_max, subavStat_dr_min, subavStat_dr_max,
        subavStat_dtP50_log10_nbins, subavStat_dr_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm2D<double> subavStat_histdttETSdr(subavStat_dttP50_log10_min, subavStat_dttP50_log10_max, subavStat_dr_min, subavStat_dr_max,
        subavStat_dttP50_log10_nbins, subavStat_dr_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm<double> subavStat_histdtP50_insidedr(subavStat_dtP50_log10_min, subavStat_dtP50_log10_max, subavStat_dtP50_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm<double> subavStat_histdtETS_insidedr(subavStat_dtP50_log10_min, subavStat_dtP50_log10_max, subavStat_dtP50_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm<double> subavStat_histdttP50_insidedr(subavStat_dttP50_log10_min, subavStat_dttP50_log10_max, subavStat_dttP50_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm<double> subavStat_histdttETS_insidedr(subavStat_dttP50_log10_min, subavStat_dttP50_log10_max, subavStat_dttP50_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm<double> subavStat_histdr_insidedr(subavStat_dr_min, subavStat_dr_max, subavStat_dr_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm<double> subavStat_histlogdr_insidedr(subavStat_dr_log10_min, subavStat_dr_log10_max, subavStat_dr_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm<double> subavStat_histdtP50_insidevtP50(subavStat_dtP50_log10_min, subavStat_dtP50_log10_max, subavStat_dtP50_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm<double> subavStat_histdtETS_insidevtP50(subavStat_dtP50_log10_min, subavStat_dtP50_log10_max, subavStat_dtP50_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm<double> subavStat_histdttP50_insidevtP50(subavStat_dttP50_log10_min, subavStat_dttP50_log10_max, subavStat_dttP50_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm<double> subavStat_histdttETS_insidevtP50(subavStat_dttP50_log10_min, subavStat_dttP50_log10_max, subavStat_dttP50_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm<double> subavStat_histdr_insidevtP50(subavStat_dr_min, subavStat_dr_max, subavStat_dr_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm<double> subavStat_histlogdr_insidevtP50(subavStat_dr_log10_min, subavStat_dr_log10_max, subavStat_dr_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm<double> subavStat_histdtP50_insidevtETS(subavStat_dtP50_log10_min, subavStat_dtP50_log10_max, subavStat_dtP50_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm<double> subavStat_histdtETS_insidevtETS(subavStat_dtP50_log10_min, subavStat_dtP50_log10_max, subavStat_dtP50_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm<double> subavStat_histdttP50_insidevtETS(subavStat_dttP50_log10_min, subavStat_dttP50_log10_max, subavStat_dttP50_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm<double> subavStat_histdttETS_insidevtETS(subavStat_dttP50_log10_min, subavStat_dttP50_log10_max, subavStat_dttP50_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm<double> subavStat_histdr_insidevtETS(subavStat_dr_min, subavStat_dr_max, subavStat_dr_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm<double> subavStat_histlogdr_insidevtETS(subavStat_dr_log10_min, subavStat_dr_log10_max, subavStat_dr_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm<double> subavStat_histdtP50_largedt(subavStat_dtP50_log10_min, subavStat_dtP50_log10_max, subavStat_dtP50_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm<double> subavStat_histdtETS_largedt(subavStat_dtP50_log10_min, subavStat_dtP50_log10_max, subavStat_dtP50_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm<double> subavStat_histdttP50_largedt(subavStat_dttP50_log10_min, subavStat_dttP50_log10_max, subavStat_dttP50_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm<double> subavStat_histdttETS_largedt(subavStat_dttP50_log10_min, subavStat_dttP50_log10_max, subavStat_dttP50_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm<double> subavStat_histdr_largedt(subavStat_dr_min, subavStat_dr_max, subavStat_dr_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm<double> subavStat_histlogdr_largedt(subavStat_dr_log10_min, subavStat_dr_log10_max, subavStat_dr_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm<double> subavStat_ijC_histdtP50(subavStat_ij_dtP50_log10_min, subavStat_ij_dtP50_log10_max, subavStat_ij_dtP50_log10_nbins, v3_math::HISTMODE_CUTOFF),
        subavStat_ijB_histdtP50(subavStat_ij_dtP50_log10_min, subavStat_ij_dtP50_log10_max, subavStat_ij_dtP50_log10_nbins, v3_math::HISTMODE_CUTOFF);
    v3_math::Histogramm<double> subavStat_ijC_histdthetaP50(subavStat_ij_dttP50_log10_min, subavStat_ij_dttP50_log10_max, subavStat_ij_dttP50_log10_nbins, v3_math::HISTMODE_CUTOFF),
        subavStat_ijB_histdthetaP50(subavStat_ij_dttP50_log10_min, subavStat_ij_dttP50_log10_max, subavStat_ij_dttP50_log10_nbins, v3_math::HISTMODE_CUTOFF);
    v3_math::WHistogramm<double> subavStat_ijWCount_whistdthetaP50(subavStat_ij_dttP50_log10_min, subavStat_ij_dttP50_log10_max, subavStat_ij_dttP50_log10_nbins, v3_math::HISTMODE_CUTOFF),
        subavStat_ijWeight_whistdthetaP50(subavStat_ij_dttP50_log10_min, subavStat_ij_dttP50_log10_max, subavStat_ij_dttP50_log10_nbins, v3_math::HISTMODE_CUTOFF);
    v3_math::Histogramm2D<double> subavStat_ij_histdtP50dr(subavStat_ij_dtP50_log10_min, subavStat_ij_dtP50_log10_max, subavStat_ij_dr_min, subavStat_ij_dr_max,
        subavStat_ij_dtP50_log10_nbins, subavStat_ij_dr_nbins, v3_math::HISTMODE_CUTOFF);
    v3_math::Histogramm2D<double> subavStat_ij_histdttP50dr(subavStat_ij_dttP50_log10_min, subavStat_ij_dttP50_log10_max, subavStat_ij_dr_min, subavStat_ij_dr_max,
        subavStat_ij_dttP50_log10_nbins, subavStat_ij_dr_nbins, v3_math::HISTMODE_CUTOFF);
    v3_math::WHistogramm2D<double> subavStat_ij_whistdtP50dr(subavStat_ij_dtP50_log10_min, subavStat_ij_dtP50_log10_max, subavStat_ij_dr_min, subavStat_ij_dr_max,
        subavStat_ij_dtP50_log10_nbins, subavStat_ij_dr_nbins, v3_math::HISTMODE_CUTOFF);
    v3_math::WHistogramm2D<double> subavStat_ij_whistdttP50dr(subavStat_ij_dttP50_log10_min, subavStat_ij_dttP50_log10_max, subavStat_ij_dr_min, subavStat_ij_dr_max,
        subavStat_ij_dttP50_log10_nbins, subavStat_ij_dr_nbins, v3_math::HISTMODE_CUTOFF);
    vector<v3_math::WHistogramm<double>> subavStat_ijCount_whistdthetaP50dr(subavStat_ij_dr_nbins), subavStat_ijWCount_whistdthetaP50dr(subavStat_ij_dr_nbins), 
        subavStat_ijWeight_whistdthetaP50dr(subavStat_ij_dr_nbins);
    for (int q = 0; q < subavStat_ij_dr_nbins; q++)
    {
        subavStat_ijCount_whistdthetaP50dr[q].change(subavStat_ij_dttP50_log10_min, subavStat_ij_dttP50_log10_max, subavStat_ij_dttP50_log10_nbins, v3_math::HISTMODE_CUTOFF);
        subavStat_ijWCount_whistdthetaP50dr[q].change(subavStat_ij_dttP50_log10_min, subavStat_ij_dttP50_log10_max, subavStat_ij_dttP50_log10_nbins, v3_math::HISTMODE_CUTOFF);
        subavStat_ijWeight_whistdthetaP50dr[q].change(subavStat_ij_dttP50_log10_min, subavStat_ij_dttP50_log10_max, subavStat_ij_dttP50_log10_nbins, v3_math::HISTMODE_CUTOFF);
    }
    double subavStat_ijCount_rbinsize = (subavStat_ij_dr_max - subavStat_ij_dr_min) / double(subavStat_ij_dr_nbins);

    for (int nf = 0; nf < nfiles; nf++)
    {
        cout << "\nFile " << nf + 1 << " / " << nfiles;
        cout.flush();

        // Sort subavalanches by timeP50
        vector<pair<double, int>> subavSorted_timeP50_subavID;
        for (int nsa = 0; nsa < subavalanches[nf].size(); nsa++)
            subavSorted_timeP50_subavID.push_back(pair<double, int>(subavalanches[nf][nsa].subavTimeP50, nsa));
        sort(subavSorted_timeP50_subavID.begin(), subavSorted_timeP50_subavID.end());

        for (int nsa = 0; nsa < subavSorted_timeP50_subavID.size(); nsa++)
        {
            // Check general time limit
            if (subavSorted_timeP50_subavID[nsa].first < analysis_time_min) continue;
            if (analysis_time_max > 0.00 && subavSorted_timeP50_subavID[nsa].first > analysis_time_max) continue;

            // Check mode
            int subavID = subavSorted_timeP50_subavID[nsa].second;
            if (simultaneousSubavalanches_mode == 1 && subavalanches_nsimultaneous[nf][subavID] != 1) continue;
            if (simultaneousSubavalanches_mode == 2 && subavalanches_nsimultaneous[nf][subavID] <= 1) continue;

            // Single subavalanche statistics
            if (subavalanches[nf][subavID].nswaps > 0.00) subavStat_histNswaps.add(log10(subavalanches[nf][subavID].nswaps));
            if (subavalanches[nf][subavID].subavDuration > 0.00) subavStat_histDuration.add(log10(subavalanches[nf][subavID].subavDuration));
            if (subavalanches[nf][subavID].subavDuration > 0.00) subavStat_whistDuration.add(log10(subavalanches[nf][subavID].subavDuration), subavalanches[nf][subavID].nswaps);
            if (subavalanches[nf][subavID].subavDiameter > 0.00) subavStat_histDiameter.add(log10(subavalanches[nf][subavID].subavDiameter));
            if (subavalanches[nf][subavID].subavRg > 0.00) subavStat_histRg.add(log10(subavalanches[nf][subavID].subavRg));
            if (subavalanches[nf][subavID].subavI2 > 0.00) subavStat_histI12.add(log10(subavalanches[nf][subavID].subavI1 / subavalanches[nf][subavID].subavI2));
            if (subavalanches[nf][subavID].subavI2 > 0.00) subavStat_whistI12.add(log10(subavalanches[nf][subavID].subavI1 / subavalanches[nf][subavID].subavI2), subavalanches[nf][subavID].nswaps);
            if (subavalanches[nf][subavID].nswaps > 0.00 && subavalanches[nf][subavID].subavDiameter > 0.00)
                subavStat_histNswapsDiameter.add(log10(subavalanches[nf][subavID].nswaps), log10(subavalanches[nf][subavID].subavDiameter));
            if (subavalanches[nf][subavID].nswaps > 0.00 && subavalanches[nf][subavID].subavRg > 0.00)
                subavStat_histNswapsRg.add(log10(subavalanches[nf][subavID].nswaps), log10(subavalanches[nf][subavID].subavRg));
            if (subavalanches[nf][subavID].nswaps > 0.00 && subavalanches[nf][subavID].subavDuration > 0.00)
                subavStat_histNswapsDuration.add(log10(subavalanches[nf][subavID].nswaps), log10(subavalanches[nf][subavID].subavDuration));
            
            // Subavalanche path analysis
            v3_phys::path_vartimestep path;
            for (int ns = 0; ns < subavalanches[nf][subavID].nswaps; ns++)
            {
                int swapID = subavalanches[nf][subavID].swaprecords_time_ID[ns].second;
                path.r.push_back(swaprecords[nf][swapID].bondCenter);
                path.t.push_back(swaprecords[nf][swapID].time);
            }
            subav_path_analyzer.add(path);

            // Inter-subavalanche statistics
            v3_math::planar::vector<double> r_this(subavalanches[nf][subavID].subavCenterX, subavalanches[nf][subavID].subavCenterY);
            for (int nsascan = nsa - 1; nsascan >= 0; nsascan--)
            {
                int prevsubavID = subavSorted_timeP50_subavID[nsa - 1].second;
                int prevNswaps = subavalanches[nf][prevsubavID].nswaps;
                if (prevNswaps >= avDetection_subavNswaps_threshold)
                {
                    double dtP50 = subavSorted_timeP50_subavID[nsa].first - subavSorted_timeP50_subavID[nsa - 1].first;
                    double dttP50 = dtP50 / subavSorted_timeP50_subavID[nsa].first;
                    double dtETS = subavalanches[nf][subavID].subavTimeStart - subavalanches[nf][prevsubavID].subavTimeEnd;
                    double dttETS = dtETS / subavalanches[nf][subavID].subavTimeStart;
                    v3_math::planar::vector<double> r_prev(subavalanches[nf][prevsubavID].subavCenterX, subavalanches[nf][prevsubavID].subavCenterY);
                    double dr = (r_this - r_prev).length();
                    double weight = v3_math::msqrt(double(subavalanches[nf][subavID].nswaps * subavalanches[nf][prevsubavID].nswaps));

                    subavStat_histdtP50.add(log10(dtP50));
                    subavStat_histdttP50.add(log10(dttP50));
                    subavStat_histdttP50_weighted.add(log10(dttP50), weight);
                    if (dtETS > 0.00) subavStat_histdtETS.add(log10(dtETS));
                    if (dttETS > 0.00) subavStat_histdttETS.add(log10(dttETS));
                    subavStat_histdr.add(dr);
                    subavStat_histdtP50dr.add(log10(dtP50), dr);
                    subavStat_histdttP50dr.add(log10(dttP50), dr);
                    if (dtETS > 0.00) subavStat_histdtETSdr.add(log10(dtETS), dr);
                    if (dttETS > 0.00) subavStat_histdttETSdr.add(log10(dttETS), dr);
                    subavStat_histlogdr.add(log10(dr));
                    subavStat_histdtP50logdr.add(log10(dtP50), log10(dr));
                    subavStat_histdttP50logdr.add(log10(dttP50), log10(dr));

                    // Large dt
                    //if (dtETS > subavalanches[nf][subavID].subavDuration && dtETS > subavalanches[nf][prevsubavID].subavDuration)
                    if (dtETS > subavalanche_dt_threshold)
                    {
                        subavStat_histdtP50_largedt.add(log10(dtP50));
                        subavStat_histdttP50_largedt.add(log10(dttP50));
                        subavStat_histdtETS_largedt.add(log10(dtETS));
                        subavStat_histdttETS_largedt.add(log10(dttETS));
                        subavStat_histdr_largedt.add(dr);
                        subavStat_histlogdr_largedt.add(log10(dr));
                    }

                    break;
                }
            }

            // Inside dr
            for (int nsascan = nsa - 1; nsascan >= 0; nsascan--)
            {
                int scansubavID = subavSorted_timeP50_subavID[nsascan].second;
                int scanNswaps = subavalanches[nf][scansubavID].nswaps;
                v3_math::planar::vector<double> r_scan(subavalanches[nf][scansubavID].subavCenterX, subavalanches[nf][scansubavID].subavCenterY);
                double dr_scan = (r_this - r_scan).length();
                if (dr_scan < avDetection_dr_threshold && scanNswaps >= avDetection_subavNswaps_threshold)
                {
                    double dtP50_insidedr = subavSorted_timeP50_subavID[nsa].first - subavSorted_timeP50_subavID[nsascan].first;
                    double dttP50_insidedr = dtP50_insidedr / subavSorted_timeP50_subavID[nsa].first;
                    double dtETS_insidedr = subavalanches[nf][subavID].subavTimeStart - subavalanches[nf][scansubavID].subavTimeEnd;
                    double dttETS_insidedr = dtETS_insidedr / subavalanches[nf][subavID].subavTimeStart;

                    subavStat_histdtP50_insidedr.add(log10(dtP50_insidedr));
                    subavStat_histdttP50_insidedr.add(log10(dttP50_insidedr));
                    if (dtETS_insidedr > 0.00) subavStat_histdtETS_insidedr.add(log10(dtETS_insidedr));
                    if (dttETS_insidedr > 0.00) subavStat_histdttETS_insidedr.add(log10(dttETS_insidedr));
                    subavStat_histdr_insidedr.add(dr_scan);
                    subavStat_histlogdr_insidedr.add(log10(dr_scan));

                    break;
                }
            }
            
            // Inside vtP50
            for (int nsascan = nsa - 1; nsascan >= 0; nsascan--)
            {
                int scansubavID = subavSorted_timeP50_subavID[nsascan].second;
                int scanNswaps = subavalanches[nf][scansubavID].nswaps;
                v3_math::planar::vector<double> r_scan(subavalanches[nf][scansubavID].subavCenterX, subavalanches[nf][scansubavID].subavCenterY);
                double dr_scan = (r_this - r_scan).length();
                double dtP50_scan = subavSorted_timeP50_subavID[nsa].first - subavSorted_timeP50_subavID[nsascan].first;
                if (dr_scan < dtP50_scan * subavStat_velocity && scanNswaps >= avDetection_subavNswaps_threshold)
                {
                    double dttP50_scan = dtP50_scan / subavSorted_timeP50_subavID[nsa].first;
                    double dtETS_scan = subavalanches[nf][subavID].subavTimeStart - subavalanches[nf][scansubavID].subavTimeEnd;
                    double dttETS_scan = dtETS_scan / subavalanches[nf][subavID].subavTimeStart;

                    subavStat_histdtP50_insidevtP50.add(log10(dtP50_scan));
                    subavStat_histdttP50_insidevtP50.add(log10(dttP50_scan));
                    if (dtETS_scan > 0.00) subavStat_histdtETS_insidevtP50.add(log10(dtETS_scan));
                    if (dttETS_scan > 0.00) subavStat_histdttETS_insidevtP50.add(log10(dttETS_scan));
                    subavStat_histdr_insidevtP50.add(dr_scan);
                    subavStat_histlogdr_insidevtP50.add(log10(dr_scan));

                    break;
                }
            }

            // Inside vtETS
            for (int nsascan = nsa - 1; nsascan >= 0; nsascan--)
            {
                int scansubavID = subavSorted_timeP50_subavID[nsascan].second;
                int scanNswaps = subavalanches[nf][scansubavID].nswaps;
                v3_math::planar::vector<double> r_scan(subavalanches[nf][scansubavID].subavCenterX, subavalanches[nf][scansubavID].subavCenterY);
                double dr_scan = (r_this - r_scan).length();
                double dtETS_scan = subavalanches[nf][subavID].subavTimeStart - subavalanches[nf][scansubavID].subavTimeEnd;
                if (dr_scan < dtETS_scan * subavStat_velocity && scanNswaps >= avDetection_subavNswaps_threshold)
                {
                    double dtP50_scan = subavSorted_timeP50_subavID[nsa].first - subavSorted_timeP50_subavID[nsascan].first;
                    double dttP50_scan = dtP50_scan / subavSorted_timeP50_subavID[nsa].first;
                    double dttETS_scan = dtETS_scan / subavalanches[nf][subavID].subavTimeStart;

                    subavStat_histdtP50_insidevtETS.add(log10(dtP50_scan));
                    subavStat_histdttP50_insidevtETS.add(log10(dttP50_scan));
                    if (dtETS_scan > 0.00) subavStat_histdtETS_insidevtETS.add(log10(dtETS_scan));
                    if (dttETS_scan > 0.00) subavStat_histdttETS_insidevtETS.add(log10(dttETS_scan));
                    subavStat_histdr_insidevtETS.add(dr_scan);
                    subavStat_histlogdr_insidevtETS.add(log10(dr_scan));

                    break;
                }
            }

            double tmax = v3_math::max(analysis_time_max, maxsimtime[nf]);
            double dtP50max = tmax - subavSorted_timeP50_subavID[nsa].first;
            subavStat_ijB_histdtP50.addAllBinsBefore(log10(dtP50max));
            double dthetaP50max = log10(tmax) - log10(subavSorted_timeP50_subavID[nsa].first);
            subavStat_ijB_histdthetaP50.addAllBinsBefore(log10(dthetaP50max));
            double dr = (subavStat_ij_dr_max - subavStat_ij_dr_min) / double(subavStat_ij_dr_nbins);
            //int swapP50ID = subavalanches[nf][subavID].swaprecords_time_ID[subavalanches[nf][subavID].nswaps / 2].second;
            //double boxsize = v3_math::msqrt(swaprecords[nf][swapP50ID].surfaceArea);
            //v3_math::planar::point<double> boxl(0.00, 0.00), boxr(boxsize, boxsize);
            //v3_math::planar::rectangle<double> rect(boxl, boxr);
            //v3_math::planar::point<double> c(r_this.x, r_this.y);
            for (int rbin = 0; rbin < subavStat_ij_dr_nbins; rbin++)
            {
                //double rmin = subavStat_ij_dr_min + double(rbin) * dr;
                //double rmax = rmin + dr;
                //v3_math::planar::circle<double> cmin(c, rmin), cmax(c, rmax);
                //double S = v3_math::planar::intersectarea(rect, cmax) - v3_math::planar::intersectarea(rect, cmin);
                //subavStat_ijWeight_whistdthetaP50dr[rbin].addAllBinsBefore(log10(dthetaP50max), v3_math::msqrt(double(subavalanches[nf][subavID].nswaps)));
                subavStat_ijWeight_whistdthetaP50.addAllBinsBefore(log10(dthetaP50max), subavalanches[nf][subavID].nswaps);
                subavStat_ijWeight_whistdthetaP50dr[rbin].addAllBinsBefore(log10(dthetaP50max), 1.00);
            }
            for (int nsa2 = nsa + 1; nsa2 < subavSorted_timeP50_subavID.size(); nsa2++)
            {
                // Check general time limit
                if (subavSorted_timeP50_subavID[nsa2].first < analysis_time_min) continue;
                if (analysis_time_max > 0.00 && subavSorted_timeP50_subavID[nsa2].first > analysis_time_max) continue;

                // Check mode
                int subavID2 = subavSorted_timeP50_subavID[nsa2].second;
                //if (simultaneousSubavalanches_mode == 1 && subavalanches_nsimultaneous[nf][subavID2] != 1) continue;
                //sif (simultaneousSubavalanches_mode == 2 && subavalanches_nsimultaneous[nf][subavID2] <= 1) continue;
                double dtP50 = subavSorted_timeP50_subavID[nsa2].first - subavSorted_timeP50_subavID[nsa].first;
                double dttP50 = dtP50 / subavSorted_timeP50_subavID[nsa2].first;
                double dthetaP50 = log10(subavSorted_timeP50_subavID[nsa2].first) - log10(subavSorted_timeP50_subavID[nsa].first);
                double dtETS = subavalanches[nf][subavID2].subavTimeStart - subavalanches[nf][subavID].subavTimeEnd;
                double dttETS = dtETS / subavSorted_timeP50_subavID[nsa2].first;
                double Wij = subavalanches[nf][subavID].nswaps * subavalanches[nf][subavID2].nswaps;
                double wnswaps = v3_math::msqrt(Wij);
                v3_math::planar::vector<double> r_nsa2(subavalanches[nf][subavID2].subavCenterX, subavalanches[nf][subavID2].subavCenterY);
                double dr = (r_this - r_nsa2).length();
                if (dtP50 > 0.00)
                {
                    subavStat_ijC_histdtP50.add(log10(dtP50));
                    subavStat_ij_histdtP50dr.add(log10(dtP50), dr);
                    subavStat_ij_whistdtP50dr.add(log10(dtP50), dr, wnswaps);
                    subavStat_ijWCount_whistdthetaP50.add(log10(dthetaP50), Wij);
                }
                if (dttP50 > 0.00)
                {
                    subavStat_ijC_histdthetaP50.add(log10(dthetaP50));
                    subavStat_ij_histdttP50dr.add(log10(dttP50), dr);
                    subavStat_ij_whistdttP50dr.add(log10(dttP50), dr, wnswaps);
                    int rbin = dr / subavStat_ijCount_rbinsize;
                    if (rbin >= 0 && rbin < subavStat_ij_dr_nbins)
                    {
                        subavStat_ijCount_whistdthetaP50dr[rbin].add(log10(dthetaP50), 1.00);
                        subavStat_ijWCount_whistdthetaP50dr[rbin].add(log10(dthetaP50), wnswaps);
                    }
                }
                /*if (dtETS > 0.00) subavStat_ij_histdtETSdr.add(log10(dtETS), dr);
                if (dttETS > 0.00) subavStat_ij_histdttETSdr.add(log10(dttETS), dr);*/
            }
        }
    }


    // INTRA AVALANCHE
    cout << "Intra-avalanche swaps\n";

    cout << "Intra-avalanche subavalanches\n";
    v3_math::Histogramm<double> subavStat_intraAv_histdtP50(subavStat_dtP50_log10_min, subavStat_dtP50_log10_max, subavStat_dtP50_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::WHistogramm<double> subavStat_intraAv_histdttP50(subavStat_dttP50_log10_min, subavStat_dttP50_log10_max, subavStat_dttP50_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::WHistogramm<double> subavStat_intraAv_histdttP50_jweighted(subavStat_dttP50_log10_min, subavStat_dttP50_log10_max, subavStat_dttP50_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::WHistogramm<double> subavStat_intraAv_histdttP50_ijweighted(subavStat_dttP50_log10_min, subavStat_dttP50_log10_max, subavStat_dttP50_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm<double> subavStat_intraAv_histdtETS(subavStat_dtP50_log10_min, subavStat_dtP50_log10_max, subavStat_dtP50_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm<double> subavStat_intraAv_histdttETS(subavStat_dttP50_log10_min, subavStat_dttP50_log10_max, subavStat_dttP50_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::WHistogramm<double> subavStat_intraAv_histdr(subavStat_dr_min, subavStat_dr_max, subavStat_dr_nbins, v3_math::HISTMODE_EXT);
    v3_math::WHistogramm<double> subavStat_intraAv_histdr_jweighted(subavStat_dr_min, subavStat_dr_max, subavStat_dr_nbins, v3_math::HISTMODE_EXT);
    v3_math::WHistogramm<double> subavStat_intraAv_histdr_ijweighted(subavStat_dr_min, subavStat_dr_max, subavStat_dr_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm<double> subavStat_intraAv_histlogdr(subavStat_dr_log10_min, subavStat_dr_log10_max, subavStat_dr_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm2D<double> subavStat_intraAv_histdtP50dr(subavStat_dtP50_log10_min, subavStat_dtP50_log10_max, subavStat_dr_min, subavStat_dr_max,
        subavStat_dtP50_log10_nbins, subavStat_dr_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm2D<double> subavStat_intraAv_histdttP50dr(subavStat_dttP50_log10_min, subavStat_dttP50_log10_max, subavStat_dr_min, subavStat_dr_max,
        subavStat_dttP50_log10_nbins, subavStat_dr_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm2D<double> subavStat_intraAv_histdtP50logdr(subavStat_dtP50_log10_min, subavStat_dtP50_log10_max, subavStat_dr_log10_min, subavStat_dr_log10_max,
        subavStat_dtP50_log10_nbins, subavStat_dr_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm2D<double> subavStat_intraAv_histdttP50logdr(subavStat_dttP50_log10_min, subavStat_dttP50_log10_max, subavStat_dr_log10_min, subavStat_dr_log10_max,
        subavStat_dttP50_log10_nbins, subavStat_dr_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm2D<double> subavStat_intraAv_histdtETSdr(subavStat_dtP50_log10_min, subavStat_dtP50_log10_max, subavStat_dr_min, subavStat_dr_max,
        subavStat_dtP50_log10_nbins, subavStat_dr_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm2D<double> subavStat_intraAv_histdttETSdr(subavStat_dttP50_log10_min, subavStat_dttP50_log10_max, subavStat_dr_min, subavStat_dr_max,
        subavStat_dttP50_log10_nbins, subavStat_dr_nbins, v3_math::HISTMODE_EXT);
    v3_math::WHistogramm<double> subavStat_intraAv_ijCount_whistdr(swapStat_ij_dr_min, swapStat_ij_dr_max, swapStat_ij_dr_nbins, v3_math::HISTMODE_EXT),
        subavStat_intraAv_ijWeight_whistdr(swapStat_ij_dr_min, swapStat_ij_dr_max, swapStat_ij_dr_nbins, v3_math::HISTMODE_EXT);
    v3_math::WHistogramm<double> subavStat_intraAv_ijCount_whistdr_alltime(swapStat_ij_dr_min, swapStat_ij_dr_max, swapStat_ij_dr_nbins, v3_math::HISTMODE_EXT),
        subavStat_intraAv_ijWeight_whistdr_alltime(swapStat_ij_dr_min, swapStat_ij_dr_max, swapStat_ij_dr_nbins, v3_math::HISTMODE_EXT);
    for (int nf = 0; nf < nfiles; nf++)
    {
        cout << "\nFile " << nf + 1 << " / " << nfiles;
        cout.flush();

        for (int nav = 0; nav < avalanches[nf].size(); nav++)
        {
            // Sort subavalanches by timeP50
            vector<pair<double, int>> subavSorted_timeP50_subavID;
            for (int nsa = 0; nsa < avalanches[nf][nav].subavalanches_timeP50_ID.size(); nsa++)
            {
                int subavID = avalanches[nf][nav].subavalanches_timeP50_ID[nsa].second;
                subavSorted_timeP50_subavID.push_back(pair<double, int>(subavalanches[nf][subavID].subavTimeP50, subavID));
            }
            sort(subavSorted_timeP50_subavID.begin(), subavSorted_timeP50_subavID.end());

            for (int nsa = 0; nsa < subavSorted_timeP50_subavID.size(); nsa++)
            {
                // Check general time limit
                if (subavSorted_timeP50_subavID[nsa].first < analysis_time_min) continue;
                if (analysis_time_max > 0.00 && subavSorted_timeP50_subavID[nsa].first > analysis_time_max) continue;

                // Check mode
                int subavID = subavSorted_timeP50_subavID[nsa].second;
                if (simultaneousSubavalanches_mode == 1 && subavalanches_nsimultaneous[nf][subavID] != 1) continue;
                if (simultaneousSubavalanches_mode == 2 && subavalanches_nsimultaneous[nf][subavID] <= 1) continue;

                // Inter-subavalanche statistics
                int swapID = subavalanches[nf][subavID].swaprecords_time_ID[0].second;
                v3_math::planar::point<double> boxLeft(system_posLeft_x, system_posLeft_y), boxRight(system_posLeft_x + swaprecords[nf][swapID].xsize, system_posLeft_y + swaprecords[nf][swapID].ysize);
                v3_math::planar::rectangle<double> box(boxLeft, boxRight);
                v3_math::planar::point<double> subavCenter(subavalanches[nf][subavID].subavCenterX, subavalanches[nf][subavID].subavCenterY);
                for (int rbin = 0; rbin < swapStat_ij_dr_nbins; rbin++)
                {
                    double rin = swapStat_ij_dr_min + double(rbin) * swapStat_ijCount_rbinsize;
                    double rout = swapStat_ij_dr_min + double(rbin + 1) * swapStat_ijCount_rbinsize;
                    v3_math::planar::circle<double> Cin(subavCenter, rin), Cout(subavCenter, rout);
                    double Sin = v3_math::planar::intersectarea(box, Cin);
                    double Sout = v3_math::planar::intersectarea(box, Cout);
                    double Sbin = Sout - Sin;
                    subavStat_intraAv_ijWeight_whistdr.add((rin + rout) * 0.50, Sbin*double(subavalanches[nf][subavID].nswaps));
                    subavStat_intraAv_ijWeight_whistdr_alltime.add((rin + rout) * 0.50, Sbin*double(subavalanches[nf][subavID].nswaps));
                }
                for (int nsaj = nsa + 1; nsaj < subavSorted_timeP50_subavID.size(); nsaj++)
                {
                    int nextsubavID = subavSorted_timeP50_subavID[nsaj].second;
                    v3_math::planar::vector<double> r_this(subavalanches[nf][subavID].subavCenterX, subavalanches[nf][subavID].subavCenterY),
                        r_next(subavalanches[nf][nextsubavID].subavCenterX, subavalanches[nf][nextsubavID].subavCenterY);
                    double dr = (r_this - r_next).length();
                    subavStat_intraAv_ijCount_whistdr.add(dr, subavalanches[nf][subavID].nswaps* subavalanches[nf][nextsubavID].nswaps);
                }
                for (int nsaj = 0; nsaj < subavSorted_timeP50_subavID.size(); nsaj++)
                {
                    if (nsaj == nsa) continue;
                    int nextsubavID = subavSorted_timeP50_subavID[nsaj].second;
                    v3_math::planar::vector<double> r_this(subavalanches[nf][subavID].subavCenterX, subavalanches[nf][subavID].subavCenterY),
                        r_next(subavalanches[nf][nextsubavID].subavCenterX, subavalanches[nf][nextsubavID].subavCenterY);
                    double dr = (r_this - r_next).length();
                    subavStat_intraAv_ijCount_whistdr_alltime.add(dr, subavalanches[nf][subavID].nswaps * subavalanches[nf][nextsubavID].nswaps);
                }
                if (nsa < 1) continue;
                int prevsubavID = subavSorted_timeP50_subavID[nsa - 1].second;
                double dtP50 = subavSorted_timeP50_subavID[nsa].first - subavSorted_timeP50_subavID[nsa - 1].first;
                double dttP50 = dtP50 / subavSorted_timeP50_subavID[nsa].first;
                double dtETS = subavalanches[nf][subavID].subavTimeStart - subavalanches[nf][prevsubavID].subavTimeEnd;
                double dttETS = dtETS / subavalanches[nf][subavID].subavTimeStart;
                v3_math::planar::vector<double> r_this(subavalanches[nf][subavID].subavCenterX, subavalanches[nf][subavID].subavCenterY),
                    r_prev(subavalanches[nf][prevsubavID].subavCenterX, subavalanches[nf][prevsubavID].subavCenterY);
                double dr = (r_this - r_prev).length();

                subavStat_intraAv_histdtP50.add(log10(dtP50));
                subavStat_intraAv_histdttP50.add(log10(dttP50), 1.00);
                subavStat_intraAv_histdttP50_jweighted.add(log10(dttP50), subavalanches[nf][subavID].nswaps);
                subavStat_intraAv_histdttP50_ijweighted.add(log10(dttP50), subavalanches[nf][subavID].nswaps*subavalanches[nf][prevsubavID].nswaps);
                if (dtETS > 0.00) subavStat_intraAv_histdtETS.add(log10(dtETS));
                if (dttETS > 0.00) subavStat_intraAv_histdttETS.add(log10(dttETS));
                subavStat_intraAv_histdr.add(dr, 1.00);
                subavStat_intraAv_histdr_jweighted.add(dr, subavalanches[nf][subavID].nswaps);
                subavStat_intraAv_histdr_ijweighted.add(dr, subavalanches[nf][subavID].nswaps* subavalanches[nf][prevsubavID].nswaps);
                subavStat_intraAv_histdtP50dr.add(log10(dtP50), dr);
                subavStat_intraAv_histdttP50dr.add(log10(dttP50), dr);
                if (dtETS > 0.00) subavStat_intraAv_histdtETSdr.add(log10(dtETS), dr);
                if (dttETS > 0.00) subavStat_intraAv_histdttETSdr.add(log10(dttETS), dr);
                subavStat_intraAv_histlogdr.add(log10(dr));
                subavStat_intraAv_histdtP50logdr.add(log10(dtP50), log10(dr));
                subavStat_intraAv_histdttP50logdr.add(log10(dttP50), log10(dr));
            }
        }
    }
    // Multiple subavalanche statistics
    cout << "Multiple subavalanches\n";
    v3_math::Histogramm<double> histMultiSubavalanches(-0.5, 10.5, 10, v3_math::HISTMODE_EXT);
    for (int nf = 0; nf < nfiles; nf++)
    {
        cout << "\nFile " << nf + 1 << " / " << nfiles;
        cout.flush();

        for (int nsa = 0; nsa < subavalanches[nf].size(); nsa++)
        {
            // Check general time limit
            if (subavalanches[nf][nsa].subavTimeP50 < analysis_time_min) continue;
            if (analysis_time_max > 0.00 && subavalanches[nf][nsa].subavTimeP50 > analysis_time_max) continue;

            histMultiSubavalanches.add(subavalanches_nsimultaneous[nf][nsa]);
        }
    }
    cout << "done.\n";
    /** / Collect subavalanche statistics **/

/**
* This block calculates statistics on thermal avalanches
**/
    /** $ Collect avalanche statistics **/
    cout << "\nCollect avalanche statistics...";
    // ALL AVALANCHES
    cout << "All avalanches\n";
    // Single avalanche statistics
    v3_math::Histogramm<double> avStat_Nsubavs(avStat_nsubavs_min, avStat_nsubavs_max, avStat_nsubavs_nbins, v3_math::HISTMODE_CUTOFF);
    v3_phys::path_analyser avStat_path_analyzer(av_path_setup), avStat_unitpath_analyzer(av_path_setup);
    vector<v3_math::Histogramm<double>> avStat_Nswaps_epoch(swapStat_nepoch* swapStat_nepochfrac + 1);
    vector<v3_math::Histogramm<double>> avStat_Duration_epoch(swapStat_nepoch* swapStat_nepochfrac + 1);
    vector<v3_math::Histogramm<double>> avStat_ThetaDuration_epoch(swapStat_nepoch* swapStat_nepochfrac + 1);
    vector<v3_math::Histogramm<double>> avStat_Tg_epoch(swapStat_nepoch* swapStat_nepochfrac + 1);
    vector<v3_math::Histogramm<double>> avStat_Thetag_epoch(swapStat_nepoch* swapStat_nepochfrac + 1);
    vector<v3_math::Histogramm<double>> avStat_Rg_epoch(swapStat_nepoch* swapStat_nepochfrac + 1);
    v3_math::Histogramm2D<double> avStat_histNswapsRg(avStat_nswaps_log10_min, avStat_nswaps_log10_max, avStat_Rg_log10_min, avStat_Rg_log10_max,
        avStat_nswaps_log10_nbins, avStat_Rg_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm2D<double> avStat_histNswapsThetaDuration(avStat_nswaps_log10_min, avStat_nswaps_log10_max,
        avStat_dtheta_log10_min, avStat_dtheta_log10_max,
        avStat_nswaps_log10_nbins, avStat_dtheta_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm2D<double> avStat_histRgThetaDuration(avStat_Rg_log10_min, avStat_Rg_log10_max,
        avStat_dtheta_log10_min, avStat_dtheta_log10_max,
        avStat_Rg_log10_nbins, avStat_dtheta_log10_nbins, v3_math::HISTMODE_EXT);
    for (int q = 0; q < swapStat_nepoch * swapStat_nepochfrac + 1; q++)
    {
        avStat_Nswaps_epoch[q].change(avStat_nswaps_log10_min, avStat_nswaps_log10_max, avStat_nswaps_log10_nbins, v3_math::HISTMODE_CUTOFF);
        avStat_Duration_epoch[q].change(avStat_duration_log10_min, avStat_duration_log10_max, avStat_duration_log10_nbins, v3_math::HISTMODE_CUTOFF);
        avStat_ThetaDuration_epoch[q].change(avStat_dtheta_log10_min, avStat_dtheta_log10_max, avStat_dtheta_log10_nbins, v3_math::HISTMODE_CUTOFF);
        avStat_Tg_epoch[q].change(avStat_tg_log10_min, avStat_tg_log10_max, avStat_tg_log10_nbins, v3_math::HISTMODE_CUTOFF);
        avStat_Thetag_epoch[q].change(avStat_dtheta_log10_min, avStat_dtheta_log10_max, avStat_dtheta_log10_nbins, v3_math::HISTMODE_CUTOFF);
        avStat_Rg_epoch[q].change(avStat_Rg_log10_min, avStat_Rg_log10_max, avStat_Rg_log10_nbins, v3_math::HISTMODE_CUTOFF);
    }
    // Inter-avalanche statistics
    v3_math::Histogramm<double> avStat_histdtP50(avStat_dtP50_log10_min, avStat_dtP50_log10_max, avStat_dtP50_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm<double> avStat_histdttP50(avStat_dttP50_log10_min, avStat_dttP50_log10_max, avStat_dttP50_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm<double> avStat_histdthetaP50(avStat_dthetaP50_log10_min, avStat_dthetaP50_log10_max, avStat_dthetaP50_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm<double> avStat_histdtETS(avStat_dtP50_log10_min, avStat_dtP50_log10_max, avStat_dtP50_log10_nbins, v3_math::HISTMODE_EXT);
    v3_math::Histogramm<double> avStat_histdttETS(avStat_dttP50_log10_min, avStat_dttP50_log10_max, avStat_dttP50_log10_nbins, v3_math::HISTMODE_EXT);
    vector<v3_math::Histogramm<double>> avStat_histdttP50_epoch(swapStat_nepoch* swapStat_nepochfrac + 1);
    vector<v3_math::Histogramm<double>> avStat_histdthetaP50_epoch(swapStat_nepoch* swapStat_nepochfrac + 1);
    vector<v3_math::Histogramm<double>> avStat_histdthetaP50_nooverlap_epoch(swapStat_nepoch* swapStat_nepochfrac + 1);
    vector<v3_math::Histogramm<double>> avStat_histdttETS_epoch(swapStat_nepoch* swapStat_nepochfrac + 1);
    v3_math::WHistogramm<double> avStat_ijCount_whistdthetaP50(avStat_ij_dtheta_log10_min, avStat_ij_dtheta_log10_max, avStat_ij_dtheta_log10_nbins, v3_math::HISTMODE_CUTOFF),
        avStat_ijWCount_whistdthetaP50(avStat_ij_dtheta_log10_min, avStat_ij_dtheta_log10_max, avStat_ij_dtheta_log10_nbins, v3_math::HISTMODE_CUTOFF),
        avStat_ijNoweight_whistdthetaP50(avStat_ij_dtheta_log10_min, avStat_ij_dtheta_log10_max, avStat_ij_dtheta_log10_nbins, v3_math::HISTMODE_CUTOFF),
        avStat_ijWeight_whistdthetaP50(avStat_ij_dtheta_log10_min, avStat_ij_dtheta_log10_max, avStat_ij_dtheta_log10_nbins, v3_math::HISTMODE_CUTOFF);
    for (int q = 0; q < swapStat_nepoch * swapStat_nepochfrac + 1; q++)
    {
        avStat_histdttP50_epoch[q].change(avStat_dttP50_log10_min, avStat_dttP50_log10_max, avStat_dttP50_log10_nbins, v3_math::HISTMODE_EXT);
        avStat_histdthetaP50_epoch[q].change(avStat_dthetaP50_log10_min, avStat_dthetaP50_log10_max, avStat_dthetaP50_log10_nbins, v3_math::HISTMODE_EXT);
        avStat_histdthetaP50_nooverlap_epoch[q].change(avStat_dthetaP50_log10_min, avStat_dthetaP50_log10_max, avStat_dthetaP50_log10_nbins, v3_math::HISTMODE_EXT);
        avStat_histdttETS_epoch[q].change(avStat_dttP50_log10_min, avStat_dttP50_log10_max, avStat_dttP50_log10_nbins, v3_math::HISTMODE_EXT);
    }
    for (int nf = 0; nf < nfiles; nf++)
    {
        cout << "\nFile " << nf + 1 << " / " << nfiles;
        cout.flush();

        // Sort avalanches by timeP50
        vector<pair<double, int>> avSorted_timeP50_avID;
        for (int na = 0; na < avalanches[nf].size(); na++)
        {
            if (avalanches[nf][na].avNswaps < avStat_nswaps_threshold) continue;
            avSorted_timeP50_avID.push_back(pair<double, int>(avalanches[nf][na].avTimeP50, na));
        }
        sort(avSorted_timeP50_avID.begin(), avSorted_timeP50_avID.end());

        for (int na = 0; na < avSorted_timeP50_avID.size(); na++)
        {
            // Check general time limit
            int avID = avSorted_timeP50_avID[na].second;
            if (avSorted_timeP50_avID[na].first < analysis_time_min) continue;
            if (analysis_time_max > 0.00 && avSorted_timeP50_avID[na].first > analysis_time_max) continue;

            // Single subavalanche statistics
            int nef = (log10(avalanches[nf][avID].avTimeP50) - swapStat_epochStart_log10) * double(swapStat_nepochfrac);
            if (avalanches[nf][avID].avNswaps > 0.00)
            {
                avStat_Nswaps_epoch[0].add(log10(avalanches[nf][avID].avNswaps));
                if (nef >= 0 && nef < swapStat_nepoch * swapStat_nepochfrac)
                    avStat_Nswaps_epoch[nef + 1].add(log10(avalanches[nf][avID].avNswaps));
            }
            if (avalanches[nf][avID].avNsubavs > 0.00) avStat_Nsubavs.add(log10(avalanches[nf][avID].avNsubavs));
            if (avalanches[nf][avID].avDuration > 0.00)
            {
                avStat_Duration_epoch[0].add(log10(avalanches[nf][avID].avDuration));
                if (nef >= 0 && nef < swapStat_nepoch * swapStat_nepochfrac)
                    avStat_Duration_epoch[nef + 1].add(log10(avalanches[nf][avID].avDuration));
            }
            if (avalanches[nf][avID].avThetaDuration > 0.00)
            {
                avStat_ThetaDuration_epoch[0].add(log10(avalanches[nf][avID].avThetaDuration));
                if (nef >= 0 && nef < swapStat_nepoch * swapStat_nepochfrac)
                    avStat_ThetaDuration_epoch[nef + 1].add(log10(avalanches[nf][avID].avThetaDuration));
            }
            if (avalanches[nf][avID].avTg > 0.00)
            {
                avStat_Tg_epoch[0].add(log10(avalanches[nf][avID].avTg));
                if (nef >= 0 && nef < swapStat_nepoch * swapStat_nepochfrac)
                    avStat_Tg_epoch[nef + 1].add(log10(avalanches[nf][avID].avTg));
            }
            if (avalanches[nf][avID].avThetag > 0.00)
            {
                avStat_Thetag_epoch[0].add(log10(avalanches[nf][avID].avThetag));
                if (nef >= 0 && nef < swapStat_nepoch * swapStat_nepochfrac)
                    avStat_Thetag_epoch[nef + 1].add(log10(avalanches[nf][avID].avThetag));
            }
            if (avalanches[nf][avID].avRg > 0.00)
            {
                avStat_Rg_epoch[0].add(log10(avalanches[nf][avID].avRg));
                if (nef >= 0 && nef < swapStat_nepoch * swapStat_nepochfrac)
                    avStat_Rg_epoch[nef + 1].add(log10(avalanches[nf][avID].avRg));
            }
            if (avalanches[nf][avID].avNswaps > 0.00 && avalanches[nf][avID].avRg > 0.00)
            {
                avStat_histNswapsRg.add(log10(avalanches[nf][avID].avNswaps), log10(avalanches[nf][avID].avRg));
            }
            if (avalanches[nf][avID].avNswaps > 0.00 && avalanches[nf][avID].avThetaDuration > 0.00)
            {
                avStat_histNswapsThetaDuration.add(log10(avalanches[nf][avID].avNswaps), log10(avalanches[nf][avID].avThetaDuration));
            }
            if (avalanches[nf][avID].avRg > 0.00 && avalanches[nf][avID].avThetaDuration > 0.00)
            {
                avStat_histRgThetaDuration.add(log10(avalanches[nf][avID].avRg), log10(avalanches[nf][avID].avThetaDuration));
            }

            // Avalanche path analysis
            v3_phys::path_vartimestep path, unitpath;
            for (int ns = 0; ns < avalanches[nf][avID].swaprecords_time_ID.size(); ns++)
            {
                int swapID = avalanches[nf][avID].swaprecords_time_ID[ns].second;
                path.r.push_back(swaprecords[nf][swapID].bondCenter);
                path.t.push_back(swaprecords[nf][swapID].time);
                unitpath.r.push_back(swaprecords[nf][swapID].bondCenter);
                unitpath.t.push_back(unitpath.r.size());
            }
            avStat_path_analyzer.add(path);
            avStat_unitpath_analyzer.add(unitpath);

            // Inter-avalanche statistics
            // i,i-1 statistics
            if (na > 0)
            {
                int prevavID = avSorted_timeP50_avID[na - 1].second;
                double dtP50 = avSorted_timeP50_avID[na].first - avSorted_timeP50_avID[na - 1].first;
                double dttP50 = dtP50 / avSorted_timeP50_avID[na].first;
                double dthetaP50 = log10(avSorted_timeP50_avID[na].first) - log10(avSorted_timeP50_avID[na - 1].first);
                double dtETS = avalanches[nf][avID].avTimeStart - avalanches[nf][prevavID].avTimeEnd;
                double dttETS = dtETS / avalanches[nf][avID].avTimeStart;
                double dthetaP50_nooverlap = 0.00;
                for (int nprev = na - 1; nprev >= 0; nprev--)
                {
                    int nprevID = avSorted_timeP50_avID[nprev].second;
                    if (avalanches[nf][avID].avTimeStart > avalanches[nf][nprevID].avTimeEnd)
                    {
                        dthetaP50_nooverlap = log10(avSorted_timeP50_avID[na].first) - log10(avSorted_timeP50_avID[nprev].first);
                        break;
                    }
                }
                avStat_histdtP50.add(log10(dtP50));
                avStat_histdttP50.add(log10(dttP50));
                avStat_histdthetaP50.add(log10(dthetaP50));
                if (dthetaP50_nooverlap > 0.00)
                    avStat_histdthetaP50_nooverlap_epoch[0].add(log10(dthetaP50_nooverlap));
                avStat_histdttP50_epoch[0].add(log10(dttP50));
                avStat_histdthetaP50_epoch[0].add(log10(dthetaP50));
                if (nef >= 0 && nef < swapStat_nepoch * swapStat_nepochfrac)
                {
                    avStat_histdttP50_epoch[nef + 1].add(log10(dttP50));
                    avStat_histdthetaP50_epoch[nef + 1].add(log10(dthetaP50));
                    if (dthetaP50_nooverlap > 0.00)
                        avStat_histdthetaP50_nooverlap_epoch[nef + 1].add(log10(dthetaP50_nooverlap));
                }
                if (dtETS > 0.00) avStat_histdtETS.add(log10(dtETS));
                if (dttETS > 0.00)
                {
                    avStat_histdttETS.add(log10(dttETS));
                    avStat_histdttETS_epoch[0].add(log10(dttETS));
                    if (nef >= 0 && nef < swapStat_nepoch * swapStat_nepochfrac)
                        avStat_histdttETS_epoch[nef + 1].add(log10(dttETS));
                }
            }
            // ij statistics
            double tmax = v3_math::max(analysis_time_max, maxsimtime[nf]);
            double dthetaP50max = log10(tmax) - log10(avalanches[nf][avID].avTimeP50);
            avStat_ijNoweight_whistdthetaP50.addAllBinsBefore(log10(dthetaP50max), 1.00);
            avStat_ijWeight_whistdthetaP50.addAllBinsBefore(log10(dthetaP50max), avalanches[nf][avID].avNswaps);
            for (int na2 = na + 1; na2 < avSorted_timeP50_avID.size(); na2++)
            {
                // Check general time limit
                if (avSorted_timeP50_avID[na2].first < analysis_time_min) continue;
                if (analysis_time_max > 0.00 && avSorted_timeP50_avID[na2].first > analysis_time_max) continue;

                int avID2 = avSorted_timeP50_avID[na2].second;
                double dthetaP50 = log10(avalanches[nf][avID2].avTimeP50) - log10(avalanches[nf][avID].avTimeP50);
                double Wij = avalanches[nf][avID].avNswaps * avalanches[nf][avID2].avNswaps;
                avStat_ijCount_whistdthetaP50.add(log10(dthetaP50), 1.00);
                avStat_ijWCount_whistdthetaP50.add(log10(dthetaP50), Wij);
            }
        }
    }
    // Multiple avalanche statistics
    vector<pair<int, int>> navalanches_nswaps;
    v3_math::Histogramm2D<double> histMultiAvalanches(0.0, 5.0, -0.5, 10.5, 20.0, 11, v3_math::HISTMODE_EXT);
    v3_math::Histogramm2D<double> histMultiAvalanchesN(0.0, 5.0, -0.5, 10.5, 20.0, 11, v3_math::HISTMODE_EXT);
    for (int nf = 0; nf < nfiles; nf++)
    {
        cout << "\nFile " << nf + 1 << " / " << nfiles;
        cout.flush();

        // Compose avalanche start and end times
        vector<pair<double, int>> avalancheStartEnd, nswaps;
        for (int na = 0; na < avalanches[nf].size(); na++)
        {
            if (avalanches[nf][na].avNswaps < 10) continue;
            if (avalanches[nf][na].avTimeStart != avalanches[nf][na].swaprecords_time_ID[0].first)
                cout << "\navTimeStart X";
            if (avalanches[nf][na].avTimeEnd != avalanches[nf][na].swaprecords_time_ID[avalanches[nf][na].swaprecords_time_ID.size() - 1].first)
                cout << "\navTimeEnd X";
            int nswapsstart = avalanches[nf][na].swaprecords_time_ID[0].second;
            int nswapsend = avalanches[nf][na].swaprecords_time_ID[avalanches[nf][na].swaprecords_time_ID.size() - 1].second;
            avalancheStartEnd.push_back(pair<double, int>(avalanches[nf][na].avTimeStart, 1));
            avalancheStartEnd.push_back(pair<double, int>(avalanches[nf][na].avTimeEnd, -1));
            nswaps.push_back(pair<double, int>(avalanches[nf][na].avTimeStart, nswapsstart));
            nswaps.push_back(pair<double, int>(avalanches[nf][na].avTimeEnd, nswapsend));
        }
        sort(avalancheStartEnd.begin(), avalancheStartEnd.end());
        sort(nswaps.begin(), nswaps.end());
        int nav = 0;
        for (int q = 0; q < avalancheStartEnd.size(); q++)
        {
            if (q > 0)
            {
                int dnswaps = nswaps[q].second - nswaps[q - 1].second;
                navalanches_nswaps.push_back(pair<int,int>(nav, dnswaps));
                if (dnswaps > 0) histMultiAvalanches.add(log10(dnswaps), nav);
                if (dnswaps > 0 && nav > 0) histMultiAvalanchesN.add(log10(dnswaps / nav), nav);
            }
            nav += avalancheStartEnd[q].second;
        }
    }
    cout << "done.\n";
    /** / Collect avalanche statistics **/


    /** $ Estimate snap likelihood **/
    cout << "\nEstimating snap likelihood...";
    v3_math::WHistogramm<double> swapStat_ij_whistdtheta_expectedRate(swapStat_ij_dtheta_min, swapStat_ij_dtheta_max, swapStat_ij_dtheta_nbins, v3_math::HISTMODE_EXT);
    double swapStat_ij_dtheta_binwidth = (swapStat_ij_dtheta_max - swapStat_ij_dtheta_min) / swapStat_ij_dtheta_nbins;
    for (int b = 0; b < swapStat_ij_dtheta_nbins; b++)
    {
        double binfrom = pow(10.00, swapStat_ij_dtheta_min + double(b) * swapStat_ij_dtheta_binwidth),
            binto = pow(10.00, swapStat_ij_dtheta_min + double(b + 1) * swapStat_ij_dtheta_binwidth);
        double log10dtheta = swapStat_ij_dtheta_min + (double(b) + 0.50) * swapStat_ij_dtheta_binwidth;
        int ndatapoints;
        double base, count;
        swapStat_ijNoweight_whistdtheta.get(log10dtheta, &ndatapoints, &base);
        swapStat_ijCount_whistdtheta.get(log10dtheta, &ndatapoints, &count);
        double rate = count / base / (binto - binfrom);
        if (base == 0.00) rate = 0.00;
        swapStat_ij_whistdtheta_expectedRate.add(log10dtheta, rate);
    }
    /*for (int nf = 0; nf < nfiles; nf++)
    {
        cout << "\nFile " << nf + 1 << " / " << nfiles;
        cout.flush();
        for (int nsi = 0; nsi < swaprecords[nf].size(); nsi++)
        {
            double Li = 0.00;
            for (int nsj = 0; nsj < nsi; nsj++)
            {
                double dtheta = log10(swaprecords[nf][nsi].time) - log10(swaprecords[nf][nsj].time);
                //double dr = (swaprecords[nf][nsi].bondCenter - swaprecords[nf][nsj].bondCenter).length();
                int npoints;
                double Lij;
                swapStat_ij_whistdtheta_expectedRate.get(log10(dtheta), &npoints, &Lij);
                if (npoints < 0)
                    cout << "!!! nsi="<<nsi<<", nsj="<<nsj<<", npoints="<<npoints<<", Lij="<<Lij<<" at dtheta="<<dtheta<<", log10(dtheta)="<< log10(dtheta)<<"\n";
                Li += Lij;
            }
            swaprecords[nf][nsi]._likelihood_ijdtheta = Li;
        }
    }*/

    cout << "done.\n";
    /** / Estimate snap likelihood **/


    /** $ Write swap data **/
    cout << "\nWriting swap data...";
    for (int nf = 0; nf < nfiles; nf++)
    {
        string logfname = swapfnames[nf].substr(0, swapfnames[nf].length() - 4) + "__markup.txt";
        ofstream outf(logfname.c_str());
        if (!outf)
        {
            cout << "\nCan not open file \"" << logfname << "\"";
            return 1;
        }
        outf << "Time\tSwap #\tBond center X\tBond center Y\tBond center Z\tSurface area\t"
            << "\tSubavalanche #\tSubavalanche size\tAvalanche #\tAvalanche size"
            << "\n";
        outf << std::fixed << std::setprecision(0);
        for (int ns = 0; ns < swaprecords[nf].size(); ns++)
        {
            outf << std::setprecision(4)
                << swaprecords[nf][ns].time 
                << std::setprecision(0)
                << "\t" << ns + 1
                << std::setprecision(4)
                << "\t" << swaprecords[nf][ns].bondCenter.x
                << "\t" << swaprecords[nf][ns].bondCenter.y
                << "\t" << swaprecords[nf][ns].bondCenter.z
                << "\t" << swaprecords[nf][ns].surfaceArea
                << std::setprecision(0)
                << "\t" << swap_subavalancheID[nf][ns]
                << "\t" << subavalanches[nf][swap_subavalancheID[nf][ns]].nswaps
                << "\t" << swaprecords_avalancheID[nf][ns]
                << "\t" << avalanches[nf][swaprecords_avalancheID[nf][ns]].avNswaps
                << "\n";
        }
        outf.close();
    }
    cout << " done.\n";
    /** / Write swap data  **/


    /** $ Write subavalanche data **/
    cout << "\nWriting subavalanche data...";
    for (int nf = 0; nf < nfiles; nf++)
    {
        string logfname = swapfnames[nf].substr(0, swapfnames[nf].length() - 4) + "__subavdata.txt";
        ofstream outf(logfname.c_str());
        if (!outf)
        {
            cout << "\nCan not open file \"" << logfname << "\"";
            return 1;
        }
        outf << "Subavalanche #\tNumber of swaps\tNumber of swaps (swap direction corrected)\tSubavalanche P50 time\tSubavalanche start time\tSubavalanche end time\tSubavalanche duration"
            << "\tSubavalanche mean rate\tSubavalanche diameter\tSubavalanche Rg\tSubavalanche I1\tSubavalanche I2\tSubavalanche mean density"
            << "\tSubavalanche X position\tSubavalanche Y position\tSimultaneous subavalanches count"
            << "\n";
        outf << std::fixed << std::setprecision(0);
        for (int nsa = 0; nsa < subavalanches[nf].size(); nsa++)
        {
            outf << std::setprecision(0) << nsa + 1
                << "\t" << subavalanches[nf][nsa].nswaps
                << "\t" << "--"
                << "\t" << std::setprecision(4) << subavalanches[nf][nsa].subavTimeP50
                << "\t" << subavalanches[nf][nsa].subavTimeStart
                << "\t" << subavalanches[nf][nsa].subavTimeEnd
                << "\t" << subavalanches[nf][nsa].subavDuration
                << "\t" << "--"
                << "\t" << std::setprecision(3) << subavalanches[nf][nsa].subavDiameter
                << "\t" << subavalanches[nf][nsa].subavRg
                << "\t" << subavalanches[nf][nsa].subavI1
                << "\t" << subavalanches[nf][nsa].subavI2
                << "\t" << "--"
                << "\t" << subavalanches[nf][nsa].subavCenterX << "\t" << subavalanches[nf][nsa].subavCenterY
                << "\t" << subavalanches_nsimultaneous[nf][nsa]
                << "\n";
        }
        outf.close();
    }
    cout << " done.\n";
    /** / Write subavalanche data  **/


    /** $ Write avalanche data  **/
    cout << "\nWriting avalanche data...";
    for (int nf = 0; nf < nfiles; nf++)
    {
        string logfname = swapfnames[nf].substr(0, swapfnames[nf].length() - 4) + "__avdata.txt";
        ofstream outf(logfname.c_str());
        if (!outf)
        {
            cout << "\nCan not open file \"" << logfname << "\"";
            return 1;
        }
        outf << "File #\tAvalanche ID\tNumber of subavalanches\tNumber of swaps\tNumber of swaps (direction corrected)"
            << "\tAvalanche start time\tAvalanche end time\tAvalanche avg time\tAvalanche P50 time"
            << "\tAvalanche duration\tAvalanche Tg"
            << "\tAvalanche Rg\tAvalanche Diameter"
            << "\n";
        outf << std::fixed << std::setprecision(0);
        for (int na = 0; na < avalanches[nf].size(); na++)
        {
            outf << std::setprecision(0) << nf + 1 << "\t" << avalanches[nf][na].avalancheID + 1
                << "\t" << avalanches[nf][na].avNsubavs << "\t" << avalanches[nf][na].avNswaps << "\t" << avalanches[nf][na].avNswapsDir;
            outf << std::setprecision(4) << "\t" << avalanches[nf][na].avTimeStart << "\t" << avalanches[nf][na].avTimeEnd
                << "\t" << avalanches[nf][na].avTimeAvg << "\t" << avalanches[nf][na].avTimeP50 << "\t" << avalanches[nf][na].avDuration
                << "\t" << avalanches[nf][na].avTg;
            outf << std::setprecision(3) << "\t" << avalanches[nf][na].avRg << "\t" << avalanches[nf][na].avDiameter;
            outf << "\n";
        }
        outf << "\n";
        outf.close();
    }
    cout << " done.\n";
    /** / Write avalanche data  **/


    /** $ Write all avalanche data  **/
    cout << "\nWriting all avalanche data...";
    string logf1name = outfname + "__avdata.txt";
    ofstream outf1(logf1name.c_str());
    if (!outf1)
    {
        cout << "\nCan not open file \"" << logf1name << "\"";
        return 1;
    }
    outf1 << "File #\tAvalanche ID\tNumber of subavalanches\tNumber of swaps\tNumber of swaps (direction corrected)"
        << "\tAvalanche start time\tAvalanche end time\tAvalanche avg time\tAvalanche P50 time"
        << "\tAvalanche duration\tAvalanche Tg"
        << "\tAvalanche Rg\tAvalanche Diameter"
        << "\n";
    outf1 << std::fixed << std::setprecision(0);
    for (int nf = 0; nf < nfiles; nf++)
    {
        for (int na = 0; na < avalanches[nf].size(); na++)
        {
            outf1 << std::setprecision(0) << nf + 1 << "\t" << avalanches[nf][na].avalancheID + 1
                << "\t" << avalanches[nf][na].avNsubavs << "\t" << avalanches[nf][na].avNswaps << "\t" << avalanches[nf][na].avNswapsDir;
            outf1 << std::setprecision(4) << "\t" << avalanches[nf][na].avTimeStart << "\t" << avalanches[nf][na].avTimeEnd
                << "\t" << avalanches[nf][na].avTimeAvg << "\t" << avalanches[nf][na].avTimeP50 << "\t" << avalanches[nf][na].avDuration
                << "\t" << avalanches[nf][na].avTg;
            outf1 << std::setprecision(3) << "\t" << avalanches[nf][na].avRg << "\t" << avalanches[nf][na].avDiameter;
            outf1 << "\n";
        }
        outf1 << "\n";
    }
    outf1.close();
    cout << " done.\n";
    /** / Write all avalanche data  **/


    /** $ Write all swap statistics **/
    cout << "\nWrite all swap statistics...";
    string logf2name = outfname + "_swaphist.txt";
    ofstream outf2(logf2name.c_str());
    if (!outf2)
    {
        cout << "\nCan not open file \"" << logf2name << "\"";
        return 1;
    }
    outf2 << "Swap dt across epochs\n";
    outf2 << "Epoch #\tEpoch start\tEpoch end\tdt avg\tdt P10\tdt P50\tdt P90\tnpoints\n";
    for (int q = 0; q < swapStat_nepoch * swapStat_nepochfrac; q++)
    {
        outf2 << q + 1 << "\t??\t??";
        if (swap_avgdt_epoch[q + 1].getcount() > 0)
            outf2 << "\t" << log10(swap_avgdt_epoch[q + 1].getavg());
        else outf2 << "\t--";
        outf2 << "\t" << swap_logdt_epoch[q + 1].getpercentile(0.10, true)
            << "\t" << swap_logdt_epoch[q + 1].getpercentile(0.50, true)
            << "\t" << swap_logdt_epoch[q + 1].getpercentile(0.90, true)
            << "\t" << swap_logdt_epoch[q + 1].getamount();
        outf2 << "\n";
    }
    outf2 << "\n\nDistribution of log t (Nswaps), count\n";
    outf2 << "log(t)";
    for (int q = 0; q < nswaps_logt.size(); q++) outf2 << "\t" << nswaps_logt[q];
    outf2 << "\n";
    v3_math::writeMultiple(swapStat_nswaps_logthist, outf2, false, 0);
    outf2 << "\n\nDistribution of swap dt\n";
    outf2 << "\n\nDistribution of log t (Nswaps), normalized\n";
    outf2 << "log(t)";
    for (int q = 0; q < nswaps_logt.size(); q++) outf2 << "\t" << nswaps_logt[q];
    outf2 << "\n";
    v3_math::writeMultiple(swapStat_nswaps_logthist, outf2, false, 3, 1.00);
    outf2 << "\n\nDistribution of swap dt\n";
    v3_math::writeMultiple(swap_logdt_epoch, outf2, false, 3, 1.00);
    outf2 << "\n\nDistribution of swap dtheta\n";
    v3_math::writeMultiple(swap_dtheta_epoch, outf2, false, 1, 1.00);
    outf2 << "\n\nDistribution of swap dtheta, PDF\n";
    v3_math::writeMultiple(swap_dtheta_epoch, outf2, false, 2, 1.00);
    outf2 << "\n\nDistribution of swap dtheta, custom bins\n";
    v3_math::writeMultiple(swap_cdtheta_epoch, outf2, false, 1, 1.00);
    outf2 << "\n\nDistribution of swap dtheta, PDF, custom bins\n";
    v3_math::writeMultiple(swap_cdtheta_epoch, outf2, false, 2, 1.00);
    outf2 << "\n\nDistribution of swap dtheta above dt threshold, custom bins\n";
    v3_math::writeMultiple(swap_cdthetaav_epoch, outf2, false, 1, 1.00);
    outf2 << "\n\nDistribution of swap dtheta above dt threshold, PDF, custom bins\n";
    v3_math::writeMultiple(swap_cdthetaav_epoch, outf2, false, 2, 1.00);
    outf2 << "\n\nDistribution of swap dt/t, custom bins\n";
    v3_math::writeMultiple(swap_cdtt_epoch, outf2, false, 1, 1.00);
    outf2 << "\n\nDistribution of swap dt/t, PDF, custom bins\n";
    v3_math::writeMultiple(swap_cdtt_epoch, outf2, false, 2, 1.00);
    outf2 << "\n\nDistribution of swap dt/t above dt threshold, custom bins\n";
    v3_math::writeMultiple(swap_cdttav_epoch, outf2, false, 1, 1.00);
    outf2 << "\n\nDistribution of swap dt/t above dt threshold, PDF, custom bins\n";
    v3_math::writeMultiple(swap_cdttav_epoch, outf2, false, 2, 1.00);
    outf2 << "\n\nDistribution of swap dt(i-1) VS dr(i-1), COUNTS\n";
    histdtdr.write(outf2, true, false);
    outf2 << "\n\nDistribution of swap dr, PDF\n";
    histdr.write(outf2, true, 0, false, 2, 1.00);
    outf2 << "\n\nDistribution of swap log(dr), PDF\n";
    histlogdr.write(outf2, true, 0, false, 3, 1.00);
    outf2 << "\n\nDistribution of swap dt(cluster) VS dr(cluster), COUNTS\n";
    histdtmindr.write(outf2, true, false);
    outf2 << "\n\nDistribution of swap log dt(i,j) VS dr(i,j), COUNTS\n";
    hist_ij_logdtdr.write(outf2, true, false);
    outf2 << "\n\nDistribution of swap log dt(i,j) VS log dr(i,j), COUNTS\n";
    hist_ij_logdtlogdr.write(outf2, true, false);
    outf2 << "\n\nDistribution of swap ij log10(dt)\n";
    outf2 << "log10(dt)\tUnweighted base\tUnweighted count\n";
    vector<v3_math::WHistogramm<double>*> swaphistlist(2);
    swaphistlist[0] = &swapStat_ijNoweight_whistdt;
    swaphistlist[1] = &swapStat_ijCount_whistdt;
    v3_math::writeMultiple(swaphistlist, outf2, false, false);
    outf2 << "\n\nDistribution of swap ij dr\n";
    outf2 << "dr\tUnweighted base\tUnweighted count\n";
    swaphistlist[0] = &swapStat_ijNoweight_whistdr;
    swaphistlist[1] = &swapStat_ijCount_whistdr;
    v3_math::writeMultiple(swaphistlist, outf2, false, false);
    outf2 << "\n\nDistribution of swap ij dr ALLTIME\n";
    outf2 << "dr\tUnweighted base\tUnweighted count\n";
    swaphistlist[0] = &swapStat_ijNoweight_whistdr_alltime;
    swaphistlist[1] = &swapStat_ijCount_whistdr_alltime;
    v3_math::writeMultiple(swaphistlist, outf2, false, false);
    outf2 << "\n\nDistribution of swap ij log10(dt) VS dr, BASE\n";
    v3_math::writeMultiple(swapStat_ijCount_whistdtdr, outf2, false, false);
    outf2 << "\n\nDistribution of swap ij log10(dt) VS dr, COUNT\n";
    v3_math::writeMultiple(swapStat_ijNoweight_whistdtdr, outf2, false, false);
    outf2 << "\n\nDistribution of swap ij log10(dtheta)\n";
    outf2 << "log10(dtheta)\tUnweighted base\tUnweighted count\n";
    //vector<v3_math::WHistogramm<double>*> swaphistlist(2);
    swaphistlist[0] = &swapStat_ijNoweight_whistdtheta;
    swaphistlist[1] = &swapStat_ijCount_whistdtheta;
    v3_math::writeMultiple(swaphistlist, outf2, false, false);
    outf2 << "\n\nDistribution of swap ij log10(dtheta) VS dr, BASE\n";
    v3_math::writeMultiple(swapStat_ijCount_whistdthetadr, outf2, false, false);
    outf2 << "\n\nDistribution of swap ij log10(dtheta) VS dr, COUNT\n";
    v3_math::writeMultiple(swapStat_ijNoweight_whistdthetadr, outf2, false, false);
    outf2.close();
    cout << "done.\n";
    /** / Write all swap statistics **/


    /** $ Write intra-subav swap statistics **/
    cout << "\nWrite intra-subav swap statistics...";
    string logf4name = outfname + "_intrasubavhist.txt";
    ofstream outf4(logf4name.c_str());
    if (!outf4)
    {
        cout << "\nCan not open file \"" << logf4name << "\"";
        return 1;
    }
    outf4 << "Intra-subavalanche swap dt across epochs\n";
    outf4 << "Epoch #\tEpoch start\tEpoch end\tdt avg\tdt P10\tdt P50\tdt P90\tnpoints\n";
    for (int q = 0; q < swapStat_nepoch * swapStat_nepochfrac; q++)
    {
        outf4 << q + 1 << "\t??\t??";
        outf4 << "\t--";
        outf4 << "\t" << intrasubav_swap_logdt_epoch[q + 1].getpercentile(0.10, true)
            << "\t" << intrasubav_swap_logdt_epoch[q + 1].getpercentile(0.50, true)
            << "\t" << intrasubav_swap_logdt_epoch[q + 1].getpercentile(0.90, true)
            << "\t" << intrasubav_swap_logdt_epoch[q + 1].getamount();
        outf4 << "\n";
    }
    outf4 << "\n\nDistribution of intra-subavalanche swap dt\n";
    v3_math::writeMultiple(intrasubav_swap_logdt_epoch, outf4, false, 3, 1.00);
    outf4 << "\n\nDistribution of intra-subavalanche swap dtheta\n";
    v3_math::writeMultiple(intrasubav_swap_dtheta_epoch, outf4, false, 1, 1.00);
    outf4 << "\n\nDistribution of intra-subavalanche swap dtheta, PDF\n";
    v3_math::writeMultiple(intrasubav_swap_dtheta_epoch, outf4, false, 2, 1.00);
    outf4 << "\n\nDistribution of intra-subavalanche swap dtheta, custom bins\n";
    v3_math::writeMultiple(intrasubav_swap_cdtheta_epoch, outf4, false, 1, 1.00);
    outf4 << "\n\nDistribution of intra-subavalanche swap dtheta, PDF, custom bins\n";
    v3_math::writeMultiple(intrasubav_swap_cdtheta_epoch, outf4, false, 2, 1.00);
    outf4 << "\n\nDistribution of intra-subavalanche swap dt/t, custom bins\n";
    v3_math::writeMultiple(intrasubav_swap_cdtt_epoch, outf4, false, 1, 1.00);
    outf4 << "\n\nDistribution of intra-subavalanche swap dt/t, PDF, custom bins\n";
    v3_math::writeMultiple(intrasubav_swap_cdtt_epoch, outf4, false, 2, 1.00);
    outf4 << "\n\nDistribution of intra-subavalanche swap dr, PDF\n";
    intrasubav_swap_histdr.write(outf4, true, 0, false, 2, 1.00);
    outf4 << "\n\nDistribution of intra-subavalanche swap log(dr), PDF\n";
    intrasubav_swap_histlogdr.write(outf4, true, 0, false, 3, 1.00);
    outf4 << "\n\nDistribution of intra-subavalanche swap closest neighbor dr, PDF\n";
    intrasubav_swap_histdrclosest.write(outf4, true, 0, false, 2, 1.00);
    outf4 << "\n\nDistribution of intra-subavalanche swap closest neighbor log(dr), PDF\n";
    intrasubav_swap_histlogdrclosest.write(outf4, true, 0, false, 3, 1.00);
    outf4 << "\n\nDistribution of intra-subavalanche swap ij dr\n";
    outf4 << "dr\tUnweighted base\tUnweighted count\n";
    swaphistlist[0] = &intrasubav_swap_ijNoweight_whistdr;
    swaphistlist[1] = &intrasubav_swap_ijCount_whistdr;
    v3_math::writeMultiple(swaphistlist, outf4, false, false);
    outf4 << "\n\nDistribution of intra-subavalanche swap ij dr ALLTIME\n";
    outf4 << "dr\tUnweighted base\tUnweighted count\n";
    swaphistlist[0] = &intrasubav_swap_ijNoweight_whistdr_alltime;
    swaphistlist[1] = &intrasubav_swap_ijCount_whistdr_alltime;
    v3_math::writeMultiple(swaphistlist, outf4, false, false);
    outf4 << "\n\n";
    //subav_path_analyzer.write(outf4);
    outf4.close();
    cout << "done.\n";
    /** / Write intra-subav swap statistics **/


    /** $ Write intra-av swap statistics **/
    cout << "\nWrite intra-av swap statistics...";
    string logf5name = outfname + "_intraavhist.txt";
    ofstream outf5(logf5name.c_str());
    if (!outf5)
    {
        cout << "\nCan not open file \"" << logf5name << "\"";
        return 1;
    }
    outf5 << "\n\nDistribution of intra-av swap dr, PDF\n";
    intraav_swap_histdr.write(outf5, true, 0, false, 2, 1.00);
    outf5 << "\n\nDistribution of intra-av swap log(dr), PDF\n";
    intraav_swap_histlogdr.write(outf5, true, 0, false, 3, 1.00);
    outf5 << "\n\nDistribution of intra-avalanche swap ij dr\n";
    outf5 << "dr\tUnweighted base\tUnweighted count\n";
    swaphistlist[0] = &intraav_swap_ijNoweight_whistdr;
    swaphistlist[1] = &intraav_swap_ijCount_whistdr;
    v3_math::writeMultiple(swaphistlist, outf5, false, false);
    outf5 << "\n\nDistribution of intra-avalanche swap ij dr ALLTIME\n";
    outf5 << "dr\tUnweighted base\tUnweighted count\n";
    swaphistlist[0] = &intraav_swap_ijNoweight_whistdr_alltime;
    swaphistlist[1] = &intraav_swap_ijCount_whistdr_alltime;
    v3_math::writeMultiple(swaphistlist, outf5, false, false);
    outf5 << "\n\n";
    //subav_path_analyzer.write(outf4);
    outf5.close();
    cout << "done.\n";
    /** / Write intra-av swap statistics **/


    /** $ Write single subavalanche statistics **/
    cout << "\nWrite single subavalanche statistics...";
    string logf6name = outfname + "_subav.txt";
    ofstream outf6(logf6name.c_str());
    if (!outf6)
    {
        cout << "\nCan not open file \"" << logf6name << "\"";
        return 1;
    }
    outf6 << "\n\n== Subavalanche distributions ==\n";
    outf6 << "\n\nDistribution of subavalanche number of swaps, COUNTS\n";
    subavStat_histNswaps.write(outf6, true, 2, false, 0);
    outf6 << "\n\nDistribution of subavalanche duration, COUNTS\n";
    subavStat_histDuration.write(outf6, true, 2, false, 0);
    outf6 << "\n\nDistribution of subavalanche duration, weighted, COUNTS\n";
    subavStat_whistDuration.write(outf6, true, 2, false);
    outf6 << "\n\nDistribution of subavalanche swap rate, COUNTS\n";
    subavStat_histSwapRate.write(outf6, true, 2, false, 0);
    outf6 << "\n\nDistribution of subavalanche diameter, COUNTS\n";
    subavStat_histDiameter.write(outf6, true, 2, false, 0);
    outf6 << "\n\nDistribution of subavalanche Rg, COUNTS\n";
    subavStat_histRg.write(outf6, true, 2, false, 0);
    outf6 << "\n\nDistribution of subavalanche I1/I2, COUNTS\n";
    subavStat_histI12.write(outf6, true, 2, false, 0);
    outf6 << "\n\nDistribution of subavalanche I1/I2, weighted, COUNTS\n";
    subavStat_whistI12.write(outf6, true, false, false);
    outf6 << "\n\nDistribution of subavalanche number of swaps VS diameter, COUNTS\n";
    subavStat_histNswapsDiameter.write(outf6, true, false);
    outf6 << "\n\nDistribution of subavalanche number of swaps VS Rg, COUNTS\n";
    subavStat_histNswapsRg.write(outf6, true, false);
    outf6 << "\n\nDistribution of subavalanche number of swaps VS Duration, COUNTS\n";
    subavStat_histNswapsDuration.write(outf6, true, false);
    outf6 << "\n\nDistribution of subavalanche number of swaps (swap direction corrected) VS diameter, COUNTS\n";
    subavStat_histNswapsDirDiameter.write(outf6, true, false);
    outf6 << "\n\nDistribution of subavalanche number of swaps (swap direction corrected) VS Rg, COUNTS\n";
    subavStat_histNswapsDirRg.write(outf6, true, false);
    outf6 << "\n\n== Subavalanche path analysis ==\n";
    subav_path_analyzer.write(outf6);
    outf6 << "\n\n== Fractal dimension analysis ==\n";
    outf6 << "Length\tAvg(N)\tAvg(log N)\tAvg(log S)\tAvg(log L)\tAvg(log S/L)\n";
    for (int q = 0; q < nlengths; q++)
    {
        outf6 << fractDimCubeLengths[q] << "\t";
        if (subav_FractDim_N_avg[q].getcount() > 0) outf6 << subav_FractDim_N_avg[q].getavg() << "\t";
        else outf6 << "--\t";
        if (subav_FractDim_logN_avg[q].getcount() > 0) outf6 << subav_FractDim_logN_avg[q].getavg() << "\t";
        else outf6 << "--\t";
        if (subav_FractDim_logS_avg[q].getcount() > 0) outf6 << subav_FractDim_logS_avg[q].getavg() << "\t";
        else outf6 << "--\t";
        if (subav_FractDim_logL_avg[q].getcount() > 0) outf6 << subav_FractDim_logL_avg[q].getavg() << "\t";
        else outf6 << "--\t";
        if (subav_FractDim_logStoL_avg[q].getcount() > 0) outf6 << subav_FractDim_logStoL_avg[q].getavg() << "\t";
        else outf6 << "--\t";
        outf6 << "\n";
    }
    outf6 << "\n\nDistribution of log(S/L)\n";
    for (int q = 0; q < nlengths; q++) outf6 << "\t" << fractDimCubeLengths[q];
    v3_math::writeMultiple(subav_FractDim_logStoL_hist, outf6, false, 3);
    outf6.close();
    cout << "done.\n";
    /** / Write single subavalanche statistics **/


    /** $ Write all avalanches fractal dimention **/
    cout << "\nWrite all avalanches fractal dimention...";
    string logf6bname = outfname + "_subavfractdim.txt";
    ofstream outf6b(logf6bname.c_str());
    if (!outf6b)
    {
        cout << "\nCan not open file \"" << logf6bname << "\"";
        return 1;
    }
    outf6b << "\n\n== Fractal dimension analysis ==\n";
    outf6b << "N(a)\n";
    outf6b << "File\tSubavalanche #\tNswaps\tTime\tRg";
    for (int q = 0; q < nlengths; q++)
        outf6b << "\t" << fractDimCubeLengths[q];
    outf6b << "\n";
    for (int nf = 0; nf < nfiles; nf++)
    {
        for (int ns = 0; ns < subavalanches[nf].size(); ns++)
        {
            if (subavalanches[nf][ns].subavTimeP50 < analysis_time_min) continue;
            if (analysis_time_max > 0.00 && subavalanches[nf][ns].subavTimeP50 > analysis_time_max) continue;
            if (subavalanches[nf][ns].nswaps < subav_fractDim_nswapsmin || subavalanches[nf][ns].nswaps > subav_fractDim_nswapsmax) continue;
            if (subavalanches[nf][ns].subavRg < subav_fractDim_Rgmin || subavalanches[nf][ns].subavRg > subav_fractDim_Rgmax) continue;

            outf6b << nf + 1 << "\t" << ns << "\t" << subavalanches[nf][ns].nswaps
                << "\t" << subavalanches[nf][ns].subavTimeP50 << "\t" << subavalanches[nf][ns].subavRg;
            for (int q = 0; q < nlengths; q++)
                outf6b << "\t" << subav_FractDim_N[nf][ns][q];
            outf6b << "\n";
        }
    }
    outf6b << "\n\n== Fractal dimension analysis - normalized ==\n";
    outf6b << "N(a)/max(N(a))\n";
    outf6b << "File\tSubavalanche #\tNswaps\tTime\tRg";
    for (int q = 0; q < nlengths; q++)
        outf6b << "\t" << fractDimCubeLengths[q];
    outf6b << "\n";
    outf6b << "Avg\t\t\t\t";
    for (int q = 0; q < nlengths; q++)
    {
        outf6b << "\t";
        if (subav_FractDim_NNmax_avg[q].getcount() > 0)  outf6b << subav_FractDim_NNmax_avg[q].getavg();
    }
    outf6b << "\n";
    for (int nf = 0; nf < nfiles; nf++)
    {
        for (int ns = 0; ns < subavalanches[nf].size(); ns++)
        {
            if (subavalanches[nf][ns].subavTimeP50 < analysis_time_min) continue;
            if (analysis_time_max > 0.00 && subavalanches[nf][ns].subavTimeP50 > analysis_time_max) continue;
            if (subavalanches[nf][ns].nswaps < subav_fractDim_nswapsmin || subavalanches[nf][ns].nswaps > subav_fractDim_nswapsmax) continue;
            if (subavalanches[nf][ns].subavRg < subav_fractDim_Rgmin || subavalanches[nf][ns].subavRg > subav_fractDim_Rgmax) continue;

            outf6b << nf + 1 << "\t" << ns << "\t" << subavalanches[nf][ns].nswaps
                << "\t" << subavalanches[nf][ns].subavTimeP50 << "\t" << subavalanches[nf][ns].subavRg;
            for (int q = 0; q < nlengths; q++)
                outf6b << "\t" << subav_FractDim_N[nf][ns][q] / subav_FractDim_N[nf][ns][nlengths - 1];
            outf6b << "\n";
        }
    }
    outf6b.close();
    cout << "done.\n";
    /** / Write all avalanches fractal dimention **/


    /** $ Write inter-subavalanche statistics **/
    cout << "\nWrite inter-subavalanche statistics...";
    string logf7name = outfname + "_inter-subav.txt";
    ofstream outf7(logf7name.c_str());
    if (!outf7)
    {
        cout << "\nCan not open file \"" << logf7name << "\"";
        return 1;
    }
    outf7 << "\n\nDistribution of inter-subavalanche log10(dt (P50)), COUNTS\n";
    subavStat_histdtP50.write(outf7, true, 2, false, 0);
    outf7 << "\n\nDistribution of inter-subavalanche log10(dt/t (P50)), COUNTS\n";
    subavStat_histdttP50.write(outf7, true, 2, false, 0);
    outf7 << "\n\nDistribution of inter-subavalanche log10(dt/t (P50)), weighted, COUNTS\n";
    subavStat_histdttP50_weighted.write(outf7, true, false, false);
    outf7 << "\n\nDistribution of inter-subavalanche log10(dt (end-start)), COUNTS\n";
    subavStat_histdtETS.write(outf7, true, 2, false, 0);
    outf7 << "\n\nDistribution of inter-subavalanche log10(dt/t (end-start)), COUNTS\n";
    subavStat_histdttETS.write(outf7, true, 2, false, 0);
    outf7 << "\n\nDistribution of inter-subavalanche dr, COUNTS\n";
    subavStat_histdr.write(outf7, true, 2, false, 0);
    outf7 << "\n\nDistribution of inter-subavalanche log10(dr), COUNTS\n";
    subavStat_histlogdr.write(outf7, true, 2, false, 0);
    outf7 << "\n\nDistribution of inter-subavalanche log10(dt (P50)) VS dr, COUNTS\n";
    subavStat_histdtP50dr.write(outf7, true, false);
    outf7 << "\n\nDistribution of inter-subavalanche log10(dt/t (P50)) VS dr, COUNTS\n";
    subavStat_histdttP50dr.write(outf7, true, false);
    outf7 << "\n\nDistribution of inter-subavalanche log10(dt (P50)) VS log10(dr), COUNTS\n";
    subavStat_histdtP50logdr.write(outf7, true, false);
    outf7 << "\n\nDistribution of inter-subavalanche log10(dt/t (P50)) VS log10(dr), COUNTS\n";
    subavStat_histdttP50logdr.write(outf7, true, false);
    outf7 << "\n\nDistribution of inter-subavalanche log10(dt (end-start)) VS dr, COUNTS\n";
    subavStat_histdtETSdr.write(outf7, true, false);
    outf7 << "\n\nDistribution of inter-subavalanche log10(dt/t (end-start)) VS dr, COUNTS\n";
    subavStat_histdttETSdr.write(outf7, true, false);

    outf7 << "\n\nDistribution of inter-subavalanche log10(dt (P50)) inside dr threshold, COUNTS\n";
    subavStat_histdtP50_insidedr.write(outf7, true, 2, false, 0);
    outf7 << "\n\nDistribution of inter-subavalanche log10(dt/t (P50)) inside dr threshold, COUNTS\n";
    subavStat_histdttP50_insidedr.write(outf7, true, 2, false, 0);
    outf7 << "\n\nDistribution of inter-subavalanche log10(dt (end-start)) inside dr threshold, COUNTS\n";
    subavStat_histdtETS_insidedr.write(outf7, true, 2, false, 0);
    outf7 << "\n\nDistribution of inter-subavalanche log10(dt/t (end-start)) inside dr threshold, COUNTS\n";
    subavStat_histdttETS_insidedr.write(outf7, true, 2, false, 0);
    outf7 << "\n\nDistribution of inter-subavalanche dr inside dr threshold, COUNTS\n";
    subavStat_histdr_insidedr.write(outf7, true, 2, false, 0);
    outf7 << "\n\nDistribution of inter-subavalanche log10(dr) inside dr threshold, COUNTS\n";
    subavStat_histlogdr_insidedr.write(outf7, true, 2, false, 0);

    outf7 << "\n\nDistribution of inter-subavalanche log10(dt (P50)) inside vtP50, COUNTS\n";
    subavStat_histdtP50_insidevtP50.write(outf7, true, 2, false, 0);
    outf7 << "\n\nDistribution of inter-subavalanche log10(dt/t (P50)) inside vtP50, COUNTS\n";
    subavStat_histdttP50_insidevtP50.write(outf7, true, 2, false, 0);
    outf7 << "\n\nDistribution of inter-subavalanche log10(dt (end-start)) inside vtP50, COUNTS\n";
    subavStat_histdtETS_insidevtP50.write(outf7, true, 2, false, 0);
    outf7 << "\n\nDistribution of inter-subavalanche log10(dt/t (end-start)) inside vtP50, COUNTS\n";
    subavStat_histdttETS_insidevtP50.write(outf7, true, 2, false, 0);
    outf7 << "\n\nDistribution of inter-subavalanche dr inside vtP50, COUNTS\n";
    subavStat_histdr_insidevtP50.write(outf7, true, 2, false, 0);
    outf7 << "\n\nDistribution of inter-subavalanche log10(dr) inside vtP50, COUNTS\n";
    subavStat_histlogdr_insidevtP50.write(outf7, true, 2, false, 0);

    outf7 << "\n\nDistribution of inter-subavalanche log10(dt (P50)) inside vtETS, COUNTS\n";
    subavStat_histdtP50_insidevtETS.write(outf7, true, 2, false, 0);
    outf7 << "\n\nDistribution of inter-subavalanche log10(dt/t (P50)) inside vtETS, COUNTS\n";
    subavStat_histdttP50_insidevtETS.write(outf7, true, 2, false, 0);
    outf7 << "\n\nDistribution of inter-subavalanche log10(dt (end-start)) inside vtETS, COUNTS\n";
    subavStat_histdtETS_insidevtETS.write(outf7, true, 2, false, 0);
    outf7 << "\n\nDistribution of inter-subavalanche log10(dt/t (end-start)) inside vtETS, COUNTS\n";
    subavStat_histdttETS_insidevtETS.write(outf7, true, 2, false, 0);
    outf7 << "\n\nDistribution of inter-subavalanche dr inside vtETS, COUNTS\n";
    subavStat_histdr_insidevtETS.write(outf7, true, 2, false, 0);
    outf7 << "\n\nDistribution of inter-subavalanche log10(dr) inside vtETS, COUNTS\n";
    subavStat_histlogdr_insidevtETS.write(outf7, true, 2, false, 0);

    outf7 << "\n\nDistribution of inter-subavalanche log10(dt (P50)) dtETS>duration, COUNTS\n";
    subavStat_histdtP50_largedt.write(outf7, true, 2, false, 0);
    outf7 << "\n\nDistribution of inter-subavalanche log10(dt/t (P50)) dtETS>duration, COUNTS\n";
    subavStat_histdttP50_largedt.write(outf7, true, 2, false, 0);
    outf7 << "\n\nDistribution of inter-subavalanche log10(dt (end-start)) dtETS>duration, COUNTS\n";
    subavStat_histdtETS_largedt.write(outf7, true, 2, false, 0);
    outf7 << "\n\nDistribution of inter-subavalanche log10(dt/t (end-start)) dtETS>duration, COUNTS\n";
    subavStat_histdttETS_largedt.write(outf7, true, 2, false, 0);
    outf7 << "\n\nDistribution of inter-subavalanche dr dtETS>duration, COUNTS\n";
    subavStat_histdr_largedt.write(outf7, true, 2, false, 0);
    outf7 << "\n\nDistribution of inter-subavalanche log10(dr) dtETS>duration, COUNTS\n";
    subavStat_histlogdr_largedt.write(outf7, true, 2, false, 0);

    outf7 << "\n\nDistribution of inter-subavalanche ij log10(dt (P50))\n";
    outf7 << "dtP50 bin center\tdtP50 from\tdtP50 to\tdt interval\tEvents count\tBin count\tAvg number of events\tAvg event freq\n";
    double dtP50min, dtP50step;
    v3_math::HISTMODE hmode;
    int np;
    subavStat_ijC_histdtP50.getparams(&dtP50min, &dtP50step, &hmode, &np);
    v3_math::Array<double> dtc, dtb;
    subavStat_ijC_histdtP50.gethist(&dtc);
    subavStat_ijB_histdtP50.gethist(&dtb);
    v3_math::Array<double> dtn(dtc.nElem(), 0.00), dtp(dtc.nElem(), 0.00);
    for (int q = 0; q < dtc.nElem(); q++)
    {
        double dtCenter = pow(10.00, dtP50min + dtP50step * (double(q) + 0.50));
        double dtFrom = pow(10.00, dtP50min + dtP50step * double(q));
        double dtTo = pow(10.00, dtP50min + dtP50step * double(q + 1));
        double dtInterval = dtTo - dtFrom;
        if (dtb(q) > 0.00)
        {
            dtn(q) = dtc(q) / dtb(q);
            dtp(q) = dtn(q) / dtInterval;
        }
        outf7 << dtCenter << "\t" << dtFrom << "\t" << dtTo << "\t" << dtInterval 
            << "\t" << dtc(q) << "\t" << dtb(q) << "\t" << dtn(q) << "\t" << dtp(q) << "\n";
    }
    outf7 << "\n\nDistribution of inter-subavalanche ij log10(dtheta (P50))\n";
    outf7 << "dthetaP50 bin center\tdthetaP50 from\tdthetaP50 to\tdtheta interval\tEvents count\tBin count\tAvg number of events\tAvg event freq\n";
    double dthetaP50min, dthetaP50step;
    subavStat_ijC_histdthetaP50.getparams(&dthetaP50min, &dthetaP50step, &hmode, &np);
    v3_math::Array<double> dthetac, dthetab;
    subavStat_ijC_histdthetaP50.gethist(&dthetac);
    subavStat_ijB_histdthetaP50.gethist(&dthetab);
    v3_math::Array<double> dthetan(dthetac.nElem(), 0.00), dthetap(dthetac.nElem(), 0.00);
    for (int q = 0; q < dthetac.nElem(); q++)
    {
        double dthetaCenter = pow(10.00, dthetaP50min + dthetaP50step * (double(q) + 0.50));
        double dthetaFrom = pow(10.00, dthetaP50min + dthetaP50step * double(q));
        double dthetaTo = pow(10.00, dthetaP50min + dthetaP50step * double(q + 1));
        double dthetaInterval = dthetaTo - dthetaFrom;
        if (dthetab(q) > 0.00)
        {
            dthetan(q) = dthetac(q) / dthetab(q);
            dthetap(q) = dthetan(q) / dthetaInterval;
        }
        outf7 << dthetaCenter << "\t" << dthetaFrom << "\t" << dthetaTo << "\t" << dthetaInterval 
            << "\t" << dthetac(q) << "\t" << dthetab(q) << "\t" << dthetan(q) << "\t" << dthetap(q) << "\n";
    }

    outf7 << "\n\nDistribution of inter-subavalanche ij log10(dt (P50)) VS dr, COUNTS\n";
    subavStat_ij_histdtP50dr.write(outf7, true, false);
    outf7 << "\n\nDistribution of inter-subavalanche ij log10(dt/t (P50)) VS dr, COUNTS\n";
    subavStat_ij_histdttP50dr.write(outf7, true, false);
    outf7 << "\n\nDistribution of inter-subavalanche weighted ij log10(dt (P50)) VS dr, COUNTS\n";
    subavStat_ij_whistdtP50dr.write(outf7, true);
    outf7 << "\n\nDistribution of inter-subavalanche weighted ij log10(dt/t (P50)) VS dr, COUNTS\n";
    subavStat_ij_whistdttP50dr.write(outf7, true);
    /*outf7 << "\n\nDistribution of inter-subavalanche ij log10(dt (end-start)) VS dr, COUNTS\n";
    subavStat_ij_histdtETSdr.write(outf7, true, false);
    outf7 << "\n\nDistribution of inter-subavalanche ij log10(dt/t (end-start)) VS dr, COUNTS\n";
    subavStat_ij_histdttETSdr.write(outf7, true, false);*/
    outf7 << "\n\nDistribution of inter-subavalanche ij log10(dtheta) VS dr, COUNTS\n";
    v3_math::writeMultiple(subavStat_ijCount_whistdthetaP50dr, outf7, false, false);
    outf7 << "\n\nDistribution of inter-subavalanche ij weighted log10(dtheta) VS dr, COUNTS\n";
    v3_math::writeMultiple(subavStat_ijWCount_whistdthetaP50dr, outf7, false, false);
    outf7 << "\n\nDistribution of inter-subavalanche ij weighted log10(dtheta) VS dr, WEIGHT\n";
    v3_math::writeMultiple(subavStat_ijWeight_whistdthetaP50dr, outf7, false, false);
    outf7 << "\n\nDistribution of inter-subavalanche ij log10(dtheta)\n";
    outf7 << "log10(dtheta)\tWeighted base\tWeighted count\n";
    vector<v3_math::WHistogramm<double>*> whistlistSubav(2);
    whistlistSubav[0] = &subavStat_ijWeight_whistdthetaP50;
    whistlistSubav[1] = &subavStat_ijWCount_whistdthetaP50;
    v3_math::writeMultiple(whistlistSubav, outf7, false, false);
    

    outf7.close();
    cout << "done.\n";
    /** / Write inter-subavalanche statistics **/


    /** $ Write intra-avalanche inter-subavalanche statistics **/
    cout << "\nWrite intra-avalanche inter-subavalanche statistics...";
    string logf7Bname = outfname + "_intraav-inter-subav.txt";
    ofstream outf7B(logf7Bname.c_str());
    if (!outf7B)
    {
        cout << "\nCan not open file \"" << logf7Bname << "\"";
        return 1;
    }
    outf7B << "\n\nDistribution of intra-avalanche inter-subavalanche log10(dt (P50)), COUNTS\n";
    subavStat_intraAv_histdtP50.write(outf7B, true, 2, false, 0);
    outf7B << "\n\nDistribution of intra-avalanche inter-subavalanche log10(dt/t (P50)), COUNTS\n";
    vector< v3_math::WHistogramm<double>> subavStat_intraAv_dttP50_hists;
    subavStat_intraAv_dttP50_hists.push_back(subavStat_intraAv_histdttP50);
    subavStat_intraAv_dttP50_hists.push_back(subavStat_intraAv_histdttP50_jweighted);
    subavStat_intraAv_dttP50_hists.push_back(subavStat_intraAv_histdttP50_ijweighted);
    v3_math::writeMultiple(subavStat_intraAv_dttP50_hists, outf7B, false, false);
    //subavStat_intraAv_histdttP50.write(outf7B, true, 2, false, 0);
    outf7B << "\n\nDistribution of intra-avalanche inter-subavalanche log10(dt (end-start)), COUNTS\n";
    subavStat_intraAv_histdtETS.write(outf7B, true, 2, false, 0);
    outf7B << "\n\nDistribution of intra-avalanche inter-subavalanche log10(dt/t (end-start)), COUNTS\n";
    subavStat_intraAv_histdttETS.write(outf7B, true, 2, false, 0);
    outf7B << "\n\nDistribution of intra-avalanche inter-subavalanche dr, COUNTS\n";
    vector< v3_math::WHistogramm<double>> subavStat_intraAv_dr_hists;
    subavStat_intraAv_dr_hists.push_back(subavStat_intraAv_histdr);
    subavStat_intraAv_dr_hists.push_back(subavStat_intraAv_histdr_jweighted);
    subavStat_intraAv_dr_hists.push_back(subavStat_intraAv_histdr_ijweighted);
    v3_math::writeMultiple(subavStat_intraAv_dr_hists, outf7B, false, false);
    //subavStat_intraAv_histdr.write(outf7B, true, 2, false, 0);
    outf7B << "\n\nDistribution of intra-avalanche inter-subavalanche log10(dr), COUNTS\n";
    subavStat_intraAv_histlogdr.write(outf7B, true, 2, false, 0);
    outf7B << "\n\nDistribution of intra-avalanche inter-subavalanche log10(dt (P50)) VS dr, COUNTS\n";
    subavStat_intraAv_histdtP50dr.write(outf7B, true, false);
    outf7B << "\n\nDistribution of intra-avalanche inter-subavalanche log10(dt/t (P50)) VS dr, COUNTS\n";
    subavStat_intraAv_histdttP50dr.write(outf7B, true, false);
    outf7B << "\n\nDistribution of intra-avalanche inter-subavalanche log10(dt (P50)) VS log10(dr), COUNTS\n";
    subavStat_intraAv_histdtP50logdr.write(outf7B, true, false);
    outf7B << "\n\nDistribution of intra-avalanche inter-subavalanche log10(dt/t (P50)) VS log10(dr), COUNTS\n";
    subavStat_intraAv_histdttP50logdr.write(outf7B, true, false);
    outf7B << "\n\nDistribution of intra-avalanche inter-subavalanche log10(dt (end-start)) VS dr, COUNTS\n";
    subavStat_intraAv_histdtETSdr.write(outf7B, true, false);
    outf7B << "\n\nDistribution of intra-avalanche inter-subavalanche log10(dt/t (end-start)) VS dr, COUNTS\n";
    subavStat_intraAv_histdttETSdr.write(outf7B, true, false);
    outf7B << "\n\nDistribution of intra-avalanche inter-subavalanche ij dr\n";
    outf7B << "dr\tWeighted base\tWeighted count\n";
    swaphistlist[0] = &subavStat_intraAv_ijWeight_whistdr;
    swaphistlist[1] = &subavStat_intraAv_ijCount_whistdr;
    v3_math::writeMultiple(swaphistlist, outf7B, false, false);
    outf7B << "\n\nDistribution of intra-avalanche inter-subavalanche ij dr ALLTIME\n";
    outf7B << "dr\tWeighted base\tWeighted count\n";
    swaphistlist[0] = &subavStat_intraAv_ijWeight_whistdr_alltime;
    swaphistlist[1] = &subavStat_intraAv_ijCount_whistdr_alltime;
    v3_math::writeMultiple(swaphistlist, outf7B, false, false);
    outf7B.close();
    cout << "done.\n";
    /** / Write intra-avalanche inter-subavalanche statistics **/


    /** $ Write all avalanches statistics **/
    cout << "\nWrite all avalanches statistics...";
    string logf8name = outfname + "_avhist.txt";
    ofstream outf8(logf8name.c_str());
    if (!outf8)
    {
        cout << "\nCan not open file \"" << logf8name << "\"";
        return 1;
    }
    outf8 << "\n\nDistribution of number of swaps, COUNTS\n";
    v3_math::writeMultiple(avStat_Nswaps_epoch, outf8, false, 0, 1.00);
    outf8 << "\n\nDistribution of number of swaps, PDF\n";
    v3_math::writeMultiple(avStat_Nswaps_epoch, outf8, false, 3, 1.00);
    outf8 << "\n\nDistribution of duration, COUNTS\n";
    v3_math::writeMultiple(avStat_Duration_epoch, outf8, false, 0, 1.00);
    outf8 << "\n\nDistribution of duration, PDF\n";
    v3_math::writeMultiple(avStat_Duration_epoch, outf8, false, 3, 1.00);
    outf8 << "\n\nDistribution of theta duration, COUNTS\n";
    v3_math::writeMultiple(avStat_ThetaDuration_epoch, outf8, false, 0, 1.00);
    outf8 << "\n\nDistribution of theta duration, PDF\n";
    v3_math::writeMultiple(avStat_ThetaDuration_epoch, outf8, false, 3, 1.00);
    outf8 << "\n\nDistribution of Tg, COUNTS\n";
    v3_math::writeMultiple(avStat_Tg_epoch, outf8, false, 0, 1.00);
    outf8 << "\n\nDistribution of Tg, PDF\n";
    v3_math::writeMultiple(avStat_Tg_epoch, outf8, false, 3, 1.00);
    outf8 << "\n\nDistribution of Thetag, COUNTS\n";
    v3_math::writeMultiple(avStat_Thetag_epoch, outf8, false, 0, 1.00);
    outf8 << "\n\nDistribution of Thetag, PDF\n";
    v3_math::writeMultiple(avStat_Thetag_epoch, outf8, false, 3, 1.00);
    outf8 << "\n\nDistribution of Rg, COUNT\n";
    v3_math::writeMultiple(avStat_Rg_epoch, outf8, false, 0, 1.00);
    outf8 << "\n\nDistribution of Rg, PDF\n";
    v3_math::writeMultiple(avStat_Rg_epoch, outf8, false, 3, 1.00);    
    outf8 << "\n\nDistribution of Nswaps vs Rg, COUNTS\n";
    avStat_histNswapsRg.write(outf8, true, false);
    outf8 << "\n\nDistribution of Nswaps vs ThetaDuration, COUNTS\n";
    avStat_histNswapsThetaDuration.write(outf8, true, false);
    outf8 << "\n\nDistribution of Rg vs ThetaDuration, COUNTS\n";
    avStat_histRgThetaDuration.write(outf8, true, false);
    outf8 << "\n\nDistribution of dt/t (P50), COUNT\n";
    v3_math::writeMultiple(avStat_histdttP50_epoch, outf8, false, 0, 1.00);
    outf8 << "\n\nDistribution of dt/t (P50), PDF\n";
    v3_math::writeMultiple(avStat_histdttP50_epoch, outf8, false, 3, 1.00);
    outf8 << "\n\nDistribution of dtheta (P50), COUNT\n";
    v3_math::writeMultiple(avStat_histdthetaP50_epoch, outf8, false, 0, 1.00);
    outf8 << "\n\nDistribution of dtheta (P50), PDF\n";
    v3_math::writeMultiple(avStat_histdthetaP50_epoch, outf8, false, 3, 1.00);
    outf8 << "\n\nDistribution of dtheta (P50), NO OVERLAP, COUNT\n";
    v3_math::writeMultiple(avStat_histdthetaP50_nooverlap_epoch, outf8, false, 0, 1.00);
    outf8 << "\n\nDistribution of dtheta (P50), NO OVERLAP, PDF\n";
    v3_math::writeMultiple(avStat_histdthetaP50_nooverlap_epoch, outf8, false, 3, 1.00);
    outf8 << "\n\nDistribution of dt/t (ETS), COUNT\n";
    v3_math::writeMultiple(avStat_histdttETS_epoch, outf8, false, 0, 1.00);
    outf8 << "\n\nDistribution of dt/t (ETS), PDF\n";
    v3_math::writeMultiple(avStat_histdttETS_epoch, outf8, false, 3, 1.00);
    outf8 << "\n\nDistribution of inter-avalanche ij log10(dtheta)\n";
    outf8 << "log10(dtheta)\tUnweighted base\tUnweighted count\tWeighted base\tWeighted count\n";
    vector<v3_math::WHistogramm<double>*> whistlist(4);
    whistlist[0] = &avStat_ijNoweight_whistdthetaP50;
    whistlist[1] = &avStat_ijCount_whistdthetaP50;
    whistlist[2] = &avStat_ijWeight_whistdthetaP50;
    whistlist[3] = &avStat_ijWCount_whistdthetaP50;
    v3_math::writeMultiple(whistlist, outf8, false, false);
    outf8 << "\n\n== Fractal dimension analysis ==\n";
    outf8 << "Length\tAvg(N)\tAvg(log N)\tAvg(log S)\tAvg(log L)\tAvg(log S/L)\n";
    for (int q = 0; q < nlengths; q++)
    {
        outf8 << fractDimCubeLengths[q] << "\t";
        if (av_FractDim_N_avg[q].getcount() > 0) outf8 << av_FractDim_N_avg[q].getavg() << "\t";
        else outf8 << "--\t";
        if (av_FractDim_logN_avg[q].getcount() > 0) outf8 << av_FractDim_logN_avg[q].getavg() << "\t";
        else outf8 << "--\t";
        if (av_FractDim_logS_avg[q].getcount() > 0) outf8 << av_FractDim_logS_avg[q].getavg() << "\t";
        else outf8 << "--\t";
        if (av_FractDim_logL_avg[q].getcount() > 0) outf8 << av_FractDim_logL_avg[q].getavg() << "\t";
        else outf8 << "--\t";
        if (av_FractDim_logStoL_avg[q].getcount() > 0) outf8 << av_FractDim_logStoL_avg[q].getavg() << "\t";
        else outf8 << "--\t";
        outf8 << "\n";
    }
    outf8 << "\n\nDistribution of log(S/L)\n";
    for (int q = 0; q < nlengths; q++) outf8 << "\t" << fractDimCubeLengths[q];
    v3_math::writeMultiple(av_FractDim_logStoL_hist, outf8, false, 3);
    outf8.close();
    cout << "done.\n";
    /** / Write all avalanches statistics **/


    /** $ Write all avalanches fractal dimention **/
    cout << "\nWrite all avalanches fractal dimention...";
    string logf9name = outfname + "_avfractdim.txt";
    ofstream outf9(logf9name.c_str());
    if (!outf9)
    {
        cout << "\nCan not open file \"" << logf9name << "\"";
        return 1;
    }
    outf9 << "\n\n== Fractal dimension analysis ==\n";
    outf9 << "N(a)\n";
    outf9 << "File\tAvalanche #\tNswaps\tTime\tRg";
    for (int q = 0; q < nlengths; q++)
        outf9 << "\t" << fractDimCubeLengths[q];
    outf9 << "\n";
    for (int nf = 0; nf < nfiles; nf++)
    {
        for (int nav = 0; nav < avalanches[nf].size(); nav++)
        {
            if (avalanches[nf][nav].avTimeP50 < analysis_time_min) continue;
            if (analysis_time_max > 0.00 && avalanches[nf][nav].avTimeP50 > analysis_time_max) continue;
            if (avalanches[nf][nav].avNswaps < av_fractDim_nswapsmin) continue;
            if (avalanches[nf][nav].avNswaps > av_fractDim_nswapsmax) continue;

            outf9 << nf + 1 << "\t" << nav << "\t" << avalanches[nf][nav].avNswaps 
                << "\t" << avalanches[nf][nav].avTimeP50 << "\t" << avalanches[nf][nav].avRg;
            for (int q = 0; q < nlengths; q++)
                outf9 << "\t" << av_FractDim_N[nf][nav][q];
            outf9 << "\n";
        }
    }
    outf9 << "\n\n== Fractal dimension analysis - normalized ==\n";
    outf9 << "N(a)/max(N(a))\n";
    outf9 << "File\tAvalanche #\tNswaps\tTime\tRg";
    for (int q = 0; q < nlengths; q++)
        outf9 << "\t" << fractDimCubeLengths[q];
    outf9 << "\n";
    for (int nf = 0; nf < nfiles; nf++)
    {
        for (int nav = 0; nav < avalanches[nf].size(); nav++)
        {
            if (avalanches[nf][nav].avTimeP50 < analysis_time_min) continue;
            if (analysis_time_max > 0.00 && avalanches[nf][nav].avTimeP50 > analysis_time_max) continue;
            if (avalanches[nf][nav].avNswaps < av_fractDim_nswapsmin) continue;
            if (avalanches[nf][nav].avNswaps > av_fractDim_nswapsmax) continue;

            outf9 << nf + 1 << "\t" << nav << "\t" << avalanches[nf][nav].avNswaps
                << "\t" << avalanches[nf][nav].avTimeP50 << "\t" << avalanches[nf][nav].avRg;
            for (int q = 0; q < nlengths; q++)
                outf9 << "\t" << av_FractDim_N[nf][nav][q] / av_FractDim_N[nf][nav][nlengths - 1];
            outf9 << "\n";
        }
    }
    outf9 << "\n\n== Fractal dimension analysis - subavalanches as points ==\n";
    outf9 << "N(a)\n";
    outf9 << "File\tAvalanche #\tNswaps\tNsubavalanches\tTime\tRg";
    for (int q = 0; q < nlengths; q++)
        outf9 << "\t" << fractDimCubeLengths[q];
    outf9 << "\n";
    for (int nf = 0; nf < nfiles; nf++)
    {
        for (int nav = 0; nav < avalanches[nf].size(); nav++)
        {
            if (avalanches[nf][nav].avTimeP50 < analysis_time_min) continue;
            if (analysis_time_max > 0.00 && avalanches[nf][nav].avTimeP50 > analysis_time_max) continue;
            if (avalanches[nf][nav].avNswaps < av_fractDim_nswapsmin) continue;
            if (avalanches[nf][nav].avNswaps > av_fractDim_nswapsmax) continue;

            outf9 << nf + 1 << "\t" << nav << "\t" << avalanches[nf][nav].avNswaps << "\t" << avalanches[nf][nav].avNsubavs
                << "\t" << avalanches[nf][nav].avTimeP50 << "\t" << avalanches[nf][nav].avRg;
            for (int q = 0; q < nlengths; q++)
                outf9 << "\t" << av_subavFractDim_N[nf][nav][q];
            outf9 << "\n";
        }
    }
    outf9 << "\n\n== Fractal dimension analysis - normalized ==\n";
    outf9 << "N(a)/max(N(a))\n";
    outf9 << "File\tAvalanche #\tNswaps\tNsubavalanches\tTime\tRg";
    for (int q = 0; q < nlengths; q++)
        outf9 << "\t" << fractDimCubeLengths[q];
    outf9 << "\n";
    for (int nf = 0; nf < nfiles; nf++)
    {
        for (int nav = 0; nav < avalanches[nf].size(); nav++)
        {
            if (avalanches[nf][nav].avTimeP50 < analysis_time_min) continue;
            if (analysis_time_max > 0.00 && avalanches[nf][nav].avTimeP50 > analysis_time_max) continue;
            if (avalanches[nf][nav].avNswaps < av_fractDim_nswapsmin) continue;
            if (avalanches[nf][nav].avNswaps > av_fractDim_nswapsmax) continue;

            outf9 << nf + 1 << "\t" << nav << "\t" << avalanches[nf][nav].avNswaps << "\t" << avalanches[nf][nav].avNsubavs
                << "\t" << avalanches[nf][nav].avTimeP50 << "\t" << avalanches[nf][nav].avRg;
            for (int q = 0; q < nlengths; q++)
                outf9 << "\t" << av_subavFractDim_N[nf][nav][q] / av_subavFractDim_N[nf][nav][nlengths - 1];
            outf9 << "\n";
        }
    }
    outf9.close();
    cout << "done.\n";
    /** / Write all avalanches fractal dimention **/


    /** $ Write multiple avalanche statistics **/
    cout << "\nWrite multiple avalanche statistics...";
    string logf10name = outfname + "_multiav.txt";
    ofstream outf10(logf10name.c_str());
    if (!outf10)
    {
        cout << "\nCan not open file \"" << logf10name << "\"";
        return 1;
    }
    outf10 << "Number of avalanches\tNumber of swaps\n";
    for (int q = 0; q < navalanches_nswaps.size(); q++)
    {
        outf10 << navalanches_nswaps[q].first << "\t" << navalanches_nswaps[q].second << "\n";
    }
    outf10 << "\n\nMap of Number of avalanches VS Number of swaps\n";
    histMultiAvalanches.write(outf10, true, false);
    outf10 << "\n\nMap of Number of avalanches VS Number of swaps / Navalanches\n";
    histMultiAvalanchesN.write(outf10, true, false);
    outf10.close();
    cout << "done.\n";
    /** / Write multiple avalanche statistics **/


    /** $ Write avtrj4 input files **/
    cout << "\nWriting avtrj4 input files...";
    for (int nf = 0; nf < nfiles; nf++)
    {
        string swapfname = swapfnames[nf].substr(0, swapfnames[nf].length() - 4) + "__avtrj4_swaprecords.txt";
        ofstream outf1(swapfname.c_str());
        if (!outf1)
        {
            cout << "\nCan not open file \"" << swapfname << "\"";
            return 1;
        }
        outf1 << "nswap\tnswapDirected\ttimestep\ttime\tdtime\ttimeRatio\tswapDirection"
            << "\txsize\tysize\tsurfaceArea\tmds\tmdsNormalized\tnsubavalanche\tbondCenter.x\tbondCenter.y\tbondCenter.z"
            << "\tdistanceToPrevSwap\tdistanceToPrevAssocSwap\tatom1\tatom2"
            << "\n";
        outf1 << std::fixed << std::setprecision(0);
        for (int nsa = 0; nsa < subavalanches[nf].size(); nsa++)
        {
            for (int q = 0; q < subavalanches[nf][nsa].swaprecords_time_ID.size(); q++)
            {
                int ns = subavalanches[nf][nsa].swaprecords_time_ID[q].second;
                outf1 << std::setprecision(0) << swaprecords[nf][ns].nswap
                    << "\t" << 0
                    << "\t" << swaprecords[nf][ns].timestep
                    << std::setprecision(4) << "\t" << swaprecords[nf][ns].time
                    << "\t" << 1.00
                    << "\t" << swaprecords[nf][ns].swapDirection
                    << "\t" << swaprecords[nf][ns].xsize
                    << "\t" << swaprecords[nf][ns].ysize
                    << "\t" << swaprecords[nf][ns].surfaceArea
                    << "\t" << 0.00
                    << "\t" << 0.00
                    << "\t" << subavalanches[nf][nsa].subavalancheID
                    << "\t" << swaprecords[nf][ns].bondCenter.x
                    << "\t" << swaprecords[nf][ns].bondCenter.y
                    << "\t" << swaprecords[nf][ns].bondCenter.z
                    << "\t" << 0.00
                    << "\t" << 0.00
                    << std::setprecision(0) << "\t" << swaprecords[nf][ns].atom1
                    << "\t" << swaprecords[nf][ns].atom2;
                outf1 << "\n";
            }
            
        }
        outf1 << "\n";

        /*string swapfname = swapfnames[nf].substr(0, swapfnames[nf].length() - 4) + "__avtrj4_swaprecords.txt";
        ofstream outf2(swapfname.c_str());
        if (!outf2)
        {
            cout << "\nCan not open file \"" << swapfname << "\"";
            return 1;
        }
        outf2 << "nsubav\tnswaps\tsubavTimeP50\tsubavTimeStart\tsubavTimeEnd\tsubavDuration\tsubavCenterX\tsubavCenterY"
            << "\tsubavRg\tsubavDiameter\tsubavDtP50\tsubavDttP50\tsubavDr\tnavalanche"
            << "\n";
        outf2 << std::fixed << std::setprecision(0);
        for (int na = 0; na < avalanches[nf].size(); na++)
        {
            outf2 << std::setprecision(0) << nf + 1 << "\t" << avalanches[nf][na].avalancheID + 1
                << "\t" << avalanches[nf][na].avNsubavs << "\t" << avalanches[nf][na].avNswaps << "\t" << avalanches[nf][na].avNswapsDir;
            outf2 << std::setprecision(4) << "\t" << avalanches[nf][na].avTimeStart << "\t" << avalanches[nf][na].avTimeEnd
                << "\t" << avalanches[nf][na].avTimeAvg << "\t" << avalanches[nf][na].avTimeP50 << "\t" << avalanches[nf][na].avDuration
                << "\t" << avalanches[nf][na].avTg;
            outf2 << std::setprecision(3) << "\t" << avalanches[nf][na].avRg << "\t" << avalanches[nf][na].avDiameter;
            outf2 << "\n";
        }
        outf2 << "\n";*/
    }
    cout << " done.\n";
    /** # Write avtrj4 input files **/

    /** $ Write swap data **/
    cout << "\nWriting swap data...";
    for (int nf = 0; nf < nfiles; nf++)
    {
        string logfname = swapfnames[nf].substr(0, swapfnames[nf].length() - 4) + "__markup.txt";
        ofstream outf(logfname.c_str());
        if (!outf)
        {
            cout << "\nCan not open file \"" << logfname << "\"";
            return 1;
        }
        outf << "Time\tSwap #\tBond center X\tBond center Y\tBond center Z\tSurface area"
            << "\tSubavalanche #\tSubavalanche size\tAvalanche #\tAvalanche size\t_likelihood_ijdtheta"
            << "\n";
        outf << std::fixed << std::setprecision(0);
        for (int ns = 0; ns < swaprecords[nf].size(); ns++)
        {
            outf << std::setprecision(4)
                << swaprecords[nf][ns].time
                << std::setprecision(0)
                << "\t" << ns + 1
                << std::setprecision(4)
                << "\t" << swaprecords[nf][ns].bondCenter.x
                << "\t" << swaprecords[nf][ns].bondCenter.y
                << "\t" << swaprecords[nf][ns].bondCenter.z
                << "\t" << swaprecords[nf][ns].surfaceArea
                << std::setprecision(0)
                << "\t" << swap_subavalancheID[nf][ns]
                << "\t" << subavalanches[nf][swap_subavalancheID[nf][ns]].nswaps
                << "\t" << swaprecords_avalancheID[nf][ns]
                << "\t" << avalanches[nf][swaprecords_avalancheID[nf][ns]].avNswaps
                << std::setprecision(4)
                << "\t" << swaprecords[nf][ns]._likelihood_ijdtheta
                << "\n";
        }
        outf.close();
    }
    cout << " done.\n";
    /** / Write swap data  **/


    /** $ Write subavalanche data **/
    cout << "\nWriting subavalanche data...";
    for (int nf = 0; nf < nfiles; nf++)
    {
        string logfname = swapfnames[nf].substr(0, swapfnames[nf].length() - 4) + "__subavdata.txt";
        ofstream outf(logfname.c_str());
        if (!outf)
        {
            cout << "\nCan not open file \"" << logfname << "\"";
            return 1;
        }
        outf << "Subavalanche #\tNumber of swaps\tNumber of swaps (swap direction corrected)\tSubavalanche P50 time\tSubavalanche start time\tSubavalanche end time\tSubavalanche duration"
            << "\tSubavalanche mean rate\tSubavalanche diameter\tSubavalanche Rg\tSubavalanche I1\tSubavalanche I2\tSubavalanche mean density"
            << "\tSubavalanche X position\tSubavalanche Y position\tSimultaneous subavalanches count"
            << "\n";
        outf << std::fixed << std::setprecision(0);
        for (int nsa = 0; nsa < subavalanches[nf].size(); nsa++)
        {
            outf << std::setprecision(0) << nsa + 1
                << "\t" << subavalanches[nf][nsa].nswaps
                << "\t" << "--"
                << "\t" << std::setprecision(4) << subavalanches[nf][nsa].subavTimeP50
                << "\t" << subavalanches[nf][nsa].subavTimeStart
                << "\t" << subavalanches[nf][nsa].subavTimeEnd
                << "\t" << subavalanches[nf][nsa].subavDuration
                << "\t" << "--"
                << "\t" << std::setprecision(3) << subavalanches[nf][nsa].subavDiameter
                << "\t" << subavalanches[nf][nsa].subavRg
                << "\t" << subavalanches[nf][nsa].subavI1
                << "\t" << subavalanches[nf][nsa].subavI2
                << "\t" << "--"
                << "\t" << subavalanches[nf][nsa].subavCenterX << "\t" << subavalanches[nf][nsa].subavCenterY
                << "\t" << subavalanches_nsimultaneous[nf][nsa]
                << "\n";
        }
        outf.close();
    }
    cout << " done.\n";
    /** / Write subavalanche data  **/


    /** $ Write avalanche data  **/
    cout << "\nWriting avalanche data...";
    for (int nf = 0; nf < nfiles; nf++)
    {
        string logfname = swapfnames[nf].substr(0, swapfnames[nf].length() - 4) + "__avdata.txt";
        ofstream outf(logfname.c_str());
        if (!outf)
        {
            cout << "\nCan not open file \"" << logfname << "\"";
            return 1;
        }
        outf << "File #\tAvalanche ID\tNumber of subavalanches\tNumber of swaps\tNumber of swaps (direction corrected)"
            << "\tAvalanche start time\tAvalanche end time\tAvalanche avg time\tAvalanche P50 time"
            << "\tAvalanche duration\tAvalanche Tg"
            << "\tAvalanche Rg\tAvalanche Diameter"
            << "\n";
        outf << std::fixed << std::setprecision(0);
        for (int na = 0; na < avalanches[nf].size(); na++)
        {
            outf << std::setprecision(0) << nf + 1 << "\t" << avalanches[nf][na].avalancheID + 1
                << "\t" << avalanches[nf][na].avNsubavs << "\t" << avalanches[nf][na].avNswaps << "\t" << avalanches[nf][na].avNswapsDir;
            outf << std::setprecision(4) << "\t" << avalanches[nf][na].avTimeStart << "\t" << avalanches[nf][na].avTimeEnd
                << "\t" << avalanches[nf][na].avTimeAvg << "\t" << avalanches[nf][na].avTimeP50 << "\t" << avalanches[nf][na].avDuration
                << "\t" << avalanches[nf][na].avTg;
            outf << std::setprecision(3) << "\t" << avalanches[nf][na].avRg << "\t" << avalanches[nf][na].avDiameter;
            outf << "\n";
        }
        outf << "\n";
        outf.close();
    }
    cout << " done.\n";
    /** / Write avalanche data  **/


    /** $ Write avalanche structure data  **/
    cout << "\nWriting avalanche structure data...";
    string logf11name = outfname + "_avstructure.txt";
    ofstream outf11(logf11name.c_str());
    if (!outf11)
    {
        cout << "\nCan not open file \"" << logf11name << "\"";
        return 1;
    }
    // Compose linear time absolute t absolute nswaps data
    //cout << "A\n";
    double tmax = 0.00;
    for (int nf = 0; nf < nfiles; nf++)
    {
        for (int na = 0; na < avalanches[nf].size(); na++)
        {
            if (avalanches[nf][na].avNswaps >= avStructure_nswaps_threshold)
            {
                double duration = avalanches[nf][na].avTimeEnd - avalanches[nf][na].avTimeStart;
                if (tmax < duration) tmax = duration;
            }
        }
    }
    //cout << "B\n";
    int ntsteps = tmax / avStructure_dt + 1; 
    vector<vector<int>> av_nswaps;
    for (int nf = 0; nf < nfiles; nf++)
    {
        //cout << nf << "\n";
        for (int na = 0; na < avalanches[nf].size(); na++)
        {
            //cout << "C" << na << "\n";
            if (avalanches[nf][na].avNswaps >= avStructure_nswaps_threshold)
            {
                //cout << "D" << na << "\n";
                vector<int> nswaps(ntsteps);
                int nt = 0;
                double t = 0.00;
                for (int ns = 0; ns < avalanches[nf][na].swaprecords_time_ID.size(); ns++)
                {
                    double time = avalanches[nf][na].swaprecords_time_ID[ns].first - avalanches[nf][na].avTimeStart;
                    while (t <= time)
                    {
                        nswaps[nt] = ns;
                        nt++;
                        t += avStructure_dt;
                    }
                }
                //cout << "E" << na << "\n";
                for (int q = nt; q < ntsteps; q++)
                    nswaps[q] = -1;
                //cout << "F" << na << "\n";
                av_nswaps.push_back(nswaps);
                //cout << "G" << na << "\n";
            }
        }
    }
    //cout << "H\n";
    outf11 << "Lin-lin absolute-absolute\n";
    outf11 << std::fixed << std::setprecision(0);
    outf11 << "Time";
    for (int q = 0; q < av_nswaps.size(); q++) outf11 << "\tAv " << q + 1;
    outf11 << "\n";
    for (int w = 0; w < ntsteps; w++)
    {
        outf11 << std::fixed << std::setprecision(3) << double(w) * avStructure_dt;
        outf11 << std::fixed << std::setprecision(0);
        for (int q = 0; q < av_nswaps.size(); q++)
        {
            if (av_nswaps[q][w] >= 0) outf11 << "\t" << av_nswaps[q][w];
            else outf11 << "\t";
        }
        outf11 << "\n";
    }
    // Compose log time absolute t absolute nswaps data
    //cout << "A\n";
    double logtmax = 0.00;
    for (int nf = 0; nf < nfiles; nf++)
    {
        for (int na = 0; na < avalanches[nf].size(); na++)
        {
            if (avalanches[nf][na].avNswaps < avStructure_nswaps_threshold) continue;
            if (avalanches[nf][na].avTimeEnd < analysis_time_min) continue;
            if (analysis_time_max > 0.00 && avalanches[nf][na].avTimeStart > analysis_time_max) continue;
            
            double logduration = log10(avalanches[nf][na].avTimeEnd) - log10(avalanches[nf][na].avTimeStart);
            if (logtmax < logduration) logtmax = logduration;
        }
    }
    //cout << "B\n";
    int nlogtsteps = logtmax / avStructure_dlogt + 1;
    vector<vector<int>> av_nswaps_logabs_linabs;
    vector<int> av_nswaps_logabs_linabs_maxnswaps;
    for (int nf = 0; nf < nfiles; nf++)
    {
        //cout << nf << "\n";
        for (int na = 0; na < avalanches[nf].size(); na++)
        {
            //cout << "C" << na << "\n";
            if (avalanches[nf][na].avNswaps < avStructure_nswaps_threshold) continue;
            if (avalanches[nf][na].avTimeEnd < analysis_time_min) continue;
            if (analysis_time_max > 0.00 && avalanches[nf][na].avTimeStart > analysis_time_max) continue;
            
            //cout << "D" << na << "\n";
            vector<int> nswaps(nlogtsteps);
            int nlogt = 0;
            double logt = 0.00;
            for (int ns = 0; ns < avalanches[nf][na].swaprecords_time_ID.size(); ns++)
            {
                double logtime = log10(avalanches[nf][na].swaprecords_time_ID[ns].first) - log10(avalanches[nf][na].avTimeStart);
                while (logt <= logtime)
                {
                    nswaps[nlogt] = ns;
                    nlogt++;
                    logt += avStructure_dlogt;
                }
            }
            //cout << "E" << na << "\n";
            for (int q = nlogt; q < nlogtsteps; q++)
                nswaps[q] = -1;
            //cout << "F" << na << "\n";
            av_nswaps_logabs_linabs.push_back(nswaps);
            av_nswaps_logabs_linabs_maxnswaps.push_back(avalanches[nf][na].avNswaps);
            //cout << "G" << na << "\n";
        }
    }
    //cout << "H\n";
    outf11 << "\n\nLog-lin absolute-absolute\n";
    outf11 << std::fixed << std::setprecision(0);
    outf11 << "Log time";
    for (int q = 0; q < av_nswaps_logabs_linabs.size(); q++) outf11 << "\tAv " << q + 1;
    outf11 << "\n";
    for (int w = 0; w < nlogtsteps; w++)
    {
        outf11 << std::fixed << std::setprecision(6) << double(w) * avStructure_dlogt;
        outf11 << std::fixed << std::setprecision(0);
        for (int q = 0; q < av_nswaps_logabs_linabs.size(); q++)
        {
            if (av_nswaps_logabs_linabs[q][w] >= 0) outf11 << "\t" << av_nswaps_logabs_linabs[q][w];
            else outf11 << "\t";
        }
        outf11 << "\n";
    }
    outf11 << "\n\nLog-lin absolute-normalized\n";
    outf11 << std::fixed << std::setprecision(0);
    outf11 << "Log time";
    for (int q = 0; q < av_nswaps_logabs_linabs.size(); q++) outf11 << "\tAv " << q + 1;
    outf11 << "\n";
    for (int w = 0; w < nlogtsteps; w++)
    {
        outf11 << std::fixed << std::setprecision(6) << double(w) * avStructure_dlogt;
        for (int q = 0; q < av_nswaps_logabs_linabs.size(); q++)
        {
            if(av_nswaps_logabs_linabs[q][w] >= 0) outf11 << "\t" << double(av_nswaps_logabs_linabs[q][w]) / double(av_nswaps_logabs_linabs_maxnswaps[q]);
            else outf11 << "\t";
        }
        outf11 << "\n";
    }
    // Compose log time normalized t absolute nswaps data
    //cout << "A\n";
    int nlogtnormsteps = 1.00 / avStructure_dlogtnorm + 1;
    vector<vector<int>> av_nswaps_lognorm_linabs;
    vector<int> av_nswaps_lognorm_linabs_maxnswaps;
    for (int nf = 0; nf < nfiles; nf++)
    {
        //cout << nf << "\n";
        for (int na = 0; na < avalanches[nf].size(); na++)
        {
            //cout << "C" << na << "\n";
            if (avalanches[nf][na].avNswaps < 2) continue;
            if (avalanches[nf][na].avNswaps < avStructure_nswaps_threshold) continue;
            if (avalanches[nf][na].avTimeEnd < analysis_time_min) continue;
            if (analysis_time_max > 0.00 && avalanches[nf][na].avTimeStart > analysis_time_max) continue;

            //cout << "D" << na << "\n";
            vector<int> nswaps(nlogtnormsteps);
            int nlogtn = 0;
            double logtn = 0.00;
            double logduration = log10(avalanches[nf][na].avTimeEnd) - log10(avalanches[nf][na].avTimeStart);
            for (int ns = 0; ns < avalanches[nf][na].swaprecords_time_ID.size(); ns++)
            {
                double logtime = log10(avalanches[nf][na].swaprecords_time_ID[ns].first) - log10(avalanches[nf][na].avTimeStart);
                double logtimen = logtime / logduration;
                while (logtn <= logtimen)
                {
                    nswaps[nlogtn] = ns;
                    nlogtn++;
                    logtn += avStructure_dlogtnorm;
                }
            }
            //cout << "E" << na << "\n";
            for (int q = nlogtn; q < nlogtnormsteps; q++)
                nswaps[q] = -1;
            //cout << "F" << na << "\n";
            av_nswaps_lognorm_linabs.push_back(nswaps);
            av_nswaps_lognorm_linabs_maxnswaps.push_back(avalanches[nf][na].avNswaps);
            //cout << "G" << na << "\n";
        }
    }
    //cout << "H\n";
    outf11 << "\n\nLog-lin normalized-absolute\n";
    outf11 << std::fixed << std::setprecision(0);
    outf11 << "Log time normalized";
    for (int q = 0; q < av_nswaps_lognorm_linabs.size(); q++) outf11 << "\tAv " << q + 1;
    outf11 << "\n";
    for (int w = 0; w < nlogtnormsteps; w++)
    {
        outf11 << std::fixed << std::setprecision(6) << double(w) * avStructure_dlogtnorm;
        outf11 << std::fixed << std::setprecision(0);
        for (int q = 0; q < av_nswaps_lognorm_linabs.size(); q++)
        {
            if (av_nswaps_lognorm_linabs[q][w] >= 0) outf11 << "\t" << av_nswaps_lognorm_linabs[q][w];
            else outf11 << "\t";
        }
        outf11 << "\n";
    }
    outf11 << "\n\nLog-lin normalized-normalized\n";
    outf11 << std::fixed << std::setprecision(0);
    outf11 << "Log time normalized";
    for (int q = 0; q < av_nswaps_lognorm_linabs.size(); q++) outf11 << "\tAv " << q + 1;
    outf11 << "\n";
    for (int w = 0; w < nlogtnormsteps; w++)
    {
        outf11 << std::fixed << std::setprecision(6) << double(w) * avStructure_dlogtnorm;
        for (int q = 0; q < av_nswaps_lognorm_linabs.size(); q++)
        {
            if (av_nswaps_lognorm_linabs[q][w] >= 0) outf11 << "\t" << double(av_nswaps_lognorm_linabs[q][w]) / double(av_nswaps_lognorm_linabs_maxnswaps[q]);
            else outf11 << "\t";
        }
        outf11 << "\n";
    }
    outf11.close();
    cout << " done.\n";
    /** / Write avalanche structure data  **/

    return 0;
}


