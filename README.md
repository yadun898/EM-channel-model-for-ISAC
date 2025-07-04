# EM-channel-model-for-ISAC

The simulation sets up fixed TX and RX antenna arrays and a synthetic target composed of multiple small PEC cylinders. Based on the excitation current at the TX side, the received field and the channel transfer matrix H are calculated using multiple scattering theory.

* MIMO_main.m — Main script:

  * Load configuration parameters
  * Generate TX/RX arrays
  * Create a target composed of small PEC cylinders
  * Compute MIMO channel matrix
  * Save results

* compute_MIMO_channel.m — Core channel computation:

  * Computes total MIMO channel H, including LoS H_I and scattered part H_s
  * Also returns scattered field E_S at receiver locations

* Supporting functions:

  * compute_incident_matrix: Computes direct path (LoS) contribution
  * compute_response_matrix: Computes wave coupling between arrays and scatterers
  * compute_Gamma_inv_PEC: Models the scattering response of PEC cylinders
  * compute_coupling_matrix: Captures multiple scattering interactions between scatterers


Authors:
Yadong Yang — PhD Candidate at 5G/6G Innovation Centre, University of Surrey, UK Email: yadong.yang@surrey.ac.uk
Prof. Gabriele Gradoni — Professor at University of Surrey Email: g.gradoni@surrey.ac.uk

Date: 2025-07-04
If you use this code for your work, please kindly cite or acknowledge the authors.


