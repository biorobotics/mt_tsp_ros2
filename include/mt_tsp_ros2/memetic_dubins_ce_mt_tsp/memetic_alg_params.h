#pragma once
class MemeticAlgParams {
  public:
    // Use values from paper or otherwise default values if paper did not specify
    MemeticAlgParams() : mutation_prob(0.1), // From paper
                         local_search_gd_step_size(0.01),
                         local_search_num_samples(20), // From paper
                         Tlp(2), // From paper
                         repair_step_size(0.01), 
                         elongation_tol(1e-2) {
    }

    // Use custom values
    MemeticAlgParams(double mutation_prob, 
                     double local_search_gd_step_size, 
                     int local_search_num_samples, 
                     int Tlp, 
                     double repair_step_size,
                     double elongation_tol) : mutation_prob(mutation_prob),
                                              local_search_gd_step_size(local_search_gd_step_size),
                                              local_search_num_samples(local_search_num_samples),
                                              Tlp(Tlp),
                                              repair_step_size(repair_step_size),
                                              elongation_tol(elongation_tol) {

    }

    double get_mutation_prob() const {
      return mutation_prob;
    }

    double get_local_search_gd_step_size() const {
      return local_search_gd_step_size;
    }

    int get_local_search_num_samples() const {
      return local_search_num_samples;
    }

    int get_Tlp() const {
      return Tlp;
    }

    double get_repair_step_size() const {
      return repair_step_size;
    }

    double get_elongation_tol() const {
      return elongation_tol;
    }

    const double mutation_prob;
    const double local_search_gd_step_size;
    const int local_search_num_samples;
    const int Tlp; // Perform local search every Tlp generations
    const double repair_step_size;
    const double elongation_tol;
};
