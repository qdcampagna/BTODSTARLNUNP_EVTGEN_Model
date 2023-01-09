# BTODSTARLNUNP_EVTGEN_Model
A modified version of EvtGen which includes a decay model to simulate new physics in the B->Dstlnu decay.

When using this Monte Carlo tool, please cite B. Bhattacharya, T. Browder, Q. Campagna, A. Datta, S. Dubey, L. Mukherjee, and A. Sibidanov, "A new tool to search for physics beyond the Standard Model in ${\bar B} \to D^* \ell^- {\bar\nu}$," arXiv:2203.07189[hep-ph]..

To make sure that the BTODSTARLNUNP decay model is working properly, use the script written in test_dstarlnu. Using the Makefile, make the run_dstarlnu script. Make sure that the TOP path listed in the Makefile is the same as where you have installed EvtGen. It should match the path listed as output when you run the "make install" command.

This generates a number of B->Dstlnu events without backgrounds. To generate a distribution, edit dstarlnu_np.dec to include whatever new physics coefficients you would like, and run the following command:

./run_dstarlnu -n #events -b Bdecay -u dstarlnu_np.dec -o filename.root

So, if you wanted to generate 10,000,000 B0 decays, you would use:

./run_dstarlnu -n 10000000 -b B0 -u dstarlnu_np.dec -o filename.root

If you want to simulate BBbar pairs in the Belle II environment from the Upsilon(4S) resonance, use the Makefile to make the run_dstarlnu2 script. To generate a distribution, edit BB_dstarlnu_np.dec to include whatever new physics coefficients you would like. The command follows exactly the same as above, but with run_dstarlnu2. As an example, if you wanted to generate 10,000,000 B0 decays, you would use:

./run_dstarlnu2 -n 10000000 -b B0 -u BB_dstarlnu_np.dec -o filename.root
