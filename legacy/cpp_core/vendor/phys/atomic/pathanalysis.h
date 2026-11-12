#pragma once
// Stub: only the setup POD is read from the INI in the core pipeline.
namespace v3_phys {

struct path_analyser_setup {
    int dr_log = 0;
    double dr_min = 0;
    double dr_max = 0;
    int dr_nbins = 0;
    int dt_log = 0;
    double dt_min = 0;
    double dt_max = 0;
    int dt_nbins = 0;
    int dpath_tau_taulog = 0;
    double dpath_tau_taumin = 0;
    double dpath_tau_taumax = 0;
    int dpath_tau_taunbins = 0;
    int dr_tau_taulog = 0;
    double dr_tau_taumin = 0;
    double dr_tau_taumax = 0;
    int dr_tau_taunbins = 0;
    int dr_dpath_dpathlog = 0;
    double dr_dpath_dpathmin = 0;
    double dr_dpath_dpathmax = 0;
    int dr_dpath_dpathnbins = 0;
};

}  // namespace v3_phys
