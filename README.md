# BTODSTARLNUNP_EVTGEN_Model
A modified version of EvtGen which includes a decay model to simulate new physics in the B->Dstlnu decay.

When using this Monte Carlo tool, please cite arXiv:2203.07189.

To make sure that the BTODSTARLNUNP decay model is working properly, use the script written in test_dstarlnu. Using the Makefile, make the run_dstarlnu script. This generates a number of B->Dstlnu events without backgrounds. To generate a distribution, edit dstarlnu_np.dec to include whatever new physics coefficients you would like, and run the following command:

./run_dstarlnu -n #events -b Bdecay -u dstarlnu_np.dec -o filename.root

So, if you wanted to generate 10,000,000 B0 decays, you would use:

./run_dstarlnu -n 10000000 -b B0 -u dstarlnu_np.dec -o filename.root
